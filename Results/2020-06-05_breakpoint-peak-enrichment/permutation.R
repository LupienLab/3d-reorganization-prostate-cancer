# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
suppressMessages(library("regioneR"))
suppressMessages(library("BSgenome.Hsapiens.UCSC.hg38.masked"))
source("../2020-02-19_chromoplexy/plotting-helper.R")


# ==============================================================================
# Data
# ==============================================================================
# load sample metadata
metadata <- fread("config.tsv", sep = "\t", header = TRUE)
SAMPLES <- metadata[, SampleID]
LABELS <- metadata[, Label]

# load SV breakpoints
breakpoints <- fread("../2020-02-19_chromoplexy/Graphs/sv-breakpoints.tsv")
breakpoint_pairs <- fread("../2020-02-19_chromoplexy/Graphs/sv-breakpoints.paired.tsv")
bp_gr <- toGRanges(breakpoints, genome = "hg38")

# load H3K27ac peaks
peaks <- rbindlist(lapply(
    SAMPLES,
    function(s) {
        dt <- fread(
            paste0("../../Data/Processed/2019-05-03_PCa-H3K27ac-peaks/Peaks/", s, "_peaks.filtered.narrowPeak"),
            sep = "\t",
            header = FALSE,
            col.names = c("chr", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak")
        )
        dt[, SampleID := s]
        return(dt)
    }
))
peaks_gr <- toGRanges(peaks, genome = "hg38")

# load ENCODE blacklist
blacklist <- fread(
    "../../Data/External/ENCODE_ChIP/ENCODE-blacklist.bed",
    sep = "\t",
    header = FALSE,
    col.names = c("chr", "start", "end")
)
blacklist_gr <- toGRanges(blacklist, genome = "hg38")

# ==============================================================================
# Analysis
# ==============================================================================
# get mask of regions to exclude from permutations
hg38 <- getGenomeAndMask("hg38")
hg38_canonical <- filterChromosomes(hg38$genome, organism="hg", chr.type = "canonical")
hg38_mask <- filterChromosomes(hg38$mask, organism="hg", chr.type = "canonical")
hg38_mask <- mergeRegions(hg38_mask, blacklist_gr)

# permutation test for number of overlaps
perm_test <- lapply(
    SAMPLES,
    function(s) {
        overlapPermTest(
            A = bp_gr[bp_gr$SampleID == s],
            B = peaks_gr[peaks_gr$SampleID == s],
            genome = hg38_canonical,
            ntimes = 1000,
            alternative = "greater",
            mask = hg38_mask
        )
    }
)
names(perm_test) <- SAMPLES

perm_test_data <- data.table(SampleID = SAMPLES)
perm_test_data[, p := apply(.SD, 1, function(r) perm_test[[r["SampleID"]]]$numOverlaps$pval)]

# local positioning adjustment to assess sensitivity to exact placement
lz <- lapply(
    SAMPLES,
    function(s) {
        localZScore(
            pt = perm_test[[s]],
            A = bp_gr[bp_gr$SampleID == s],
            B = peaks_gr[peaks_gr$SampleID == s],
            ntimes = 100,
            window = 100000,
            step = 10000
        )
    }
)
names(lz) <- SAMPLES

# sv overlaps
bp_peaks <- rbindlist(lapply(
    SAMPLES,
    function(s) {
        # get the sample-specific peaks and breakpoints
        sample_breaks <- bp_gr[bp_gr$SampleID == s]
        sample_peaks <- peaks_gr[peaks_gr$SampleID == s]
        # find overlaps between them
        overlaps <- as.data.table(findOverlaps(sample_breaks, sample_peaks))
        overlaps[, SampleID := s]
        # assign the breakpoint_ID to the relevant queryHit ID
        overlaps[, breakpoint_ID := apply(.SD, 1, function(r) {sample_breaks[as.numeric(r["queryHits"])]$breakpoint_ID})]
        return(overlaps)
    }
))

# summarize number of peaks overlapped for each SV breakpoint
bp_peaks_counts <- bp_peaks[, .N, by = c("SampleID", "breakpoint_ID")]

# add breakpoints that did not overlap any peaks
## start by getting the breakpoint IDs that overlapped 0 peaks
no_peaks <- data.table(
    breakpoint_ID = setdiff(breakpoints[, breakpoint_ID], bp_peaks_counts[, breakpoint_ID]),
    N = 0
)
## merge SampleIDs from the breakpoint IDs
no_peaks <- merge(
    x = no_peaks[, .SD, keyby = "breakpoint_ID"],
    y = breakpoints[, .SD, .SDcols = c("SampleID", "breakpoint_ID")],
    by = "breakpoint_ID"
)[, .SD, .SDcols = c("SampleID", "breakpoint_ID", "N")]
## concatenate results
bp_peaks_counts <- rbindlist(list(bp_peaks_counts, no_peaks))

# merge metadata information for plotting
bp_peaks_counts <- merge(
    x = bp_peaks_counts,
    y = metadata,
    by = "SampleID"
)

# ==============================================================================
# Plots
# ==============================================================================
for (s in SAMPLES) {
    cat(s, "\n")
    png(paste0("Plots/", s, ".permutation.png"), width = 20, height = 12, units = "cm", res = 300)
    plot(perm_test[[s]])
    dev.off()
    
    png(paste0("Plots/", s, ".local-z.png"), width = 20, height = 12, units = "cm", res = 300)
    plot(lz[[s]])
    dev.off()
}

# counts of how many H3K27ac peaks each SV breakpoint overlaps
gg <- (
    ggplot(data = bp_peaks_counts)
    + geom_col(
        aes(
            x = factor(breakpoint_ID, ordered = TRUE, levels = bp_peaks_counts[order(Label, N), breakpoint_ID]),
            y = N,
            group = Label,
            fill = Label
        ),
        position = "dodge"
    )
    + labs(x = "Breakpoints", y = "H3K27ac Peaks Overlapping Each Breakpoint")
    + scale_fill_manual(
        breaks = metadata[, Label],
        labels = metadata[, Label],
        values = metadata[, Sample_Colour],
        name = "Patient"
    )
    + theme_minimal()
    + theme(
        legend.position = "bottom",
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank()
    )
)
savefig(gg, "Plots/bp_peaks_counts", width = 30)

# p-value histogram for permutation tests
gg_p <- (
    ggplot(data = perm_test_data)
    + geom_col(aes(x = SampleID, y = -log10(p), fill = SampleID))
    + geom_hline(aes(yintercept = -log10(0.05)), linetype = "dashed")
    + scale_x_discrete(
        breaks = metadata[, SampleID],
        labels = metadata[, Label],
        name = NULL
    )
    + scale_y_continuous(
        minor_breaks = -log10(c(seq(1, 10) * 1e-1, seq(1, 9) * 1e-2, seq(1, 9) * 1e-3)),
        breaks = -log10(c(1, 0.1, 0.05, 0.01, 0.001)),
        labels = c(1, 0.1, 0.05, 0.01, 0.001),
        name = "p-value"
    )
    + scale_fill_manual(
        breaks = metadata[, SampleID],
        labels = metadata[, Label],
        values = metadata[, Sample_Colour]
    )
    + guides(fill = FALSE)
    + theme_minimal()
    + theme(
        legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 1, hjust = 0.5)
    )
)
savefig(gg_p, "Plots/permutation.p-values")


# ==============================================================================
# Save data
# ==============================================================================
fwrite(
    bp_peaks_counts[, .SD, .SDcols = c("N"), keyby = c("SampleID", "breakpoint_ID")],
    "Overlaps/sv-breakpoints.overlaps.tsv",
    sep = "\t",
    col.names = TRUE
)

fwrite(
    perm_test_data,
    "Overlaps/sv-breakpoints.permutations.tsv",
    sep = "\t",
    col.names = TRUE
)

saveRDS(perm_test, "Overlaps/permutation-tests.rds")
saveRDS(lz, "Overlaps/local-dependency.rds")
