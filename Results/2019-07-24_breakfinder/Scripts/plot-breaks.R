# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("RCircos"))
suppressMessages(library("ggplot2"))

CHRS = paste0('chr', c(1:22, "X", "Y"))

hg38 = data.table(
    Chrom = factor(CHRS, levels = CHRS, ordered = TRUE),
    Length = c(
        248956422,
        242193529,
        198295559,
        190214555,
        181538259,
        170805979,
        159345973,
        145138636,
        138394717,
        133797422,
        135086622,
        133275309,
        114364328,
        107043718,
        101991189,
        90338345,
        83257441,
        80373285,
        58617616,
        64444167,
        46709983,
        50818468,
        156040895,
        57227415
    )
)
hg38[, Bins := ceiling(Length / 10^6)]

# ==============================================================================
# Data
# ==============================================================================
# load metadata
metadata <- fread("../../Data/External/LowC_Samples_Data_Available.tsv", sep = "\t", header = TRUE)
SAMPLES <- paste0("PCa", metadata[, get("Sample ID")])

# load breakpoint data
breakpoints <- rbindlist(lapply(
    SAMPLES,
    function(s) {
        dt <- fread(
            paste0("Breakpoints/Default/", s, ".breaks.sorted.manually-resolved.tsv"),
            sep = "\t",
            header = FALSE,
            col.names = c(
                "chr_from", "start_from", "end_from", "chr_to", "start_to", "end_to",
                "name", "score", "strand_from", "strand_to", "resolution",
                "type", "notes"
            )
        )
        dt[, Sample := s]
        return(dt)
    }
))
# remove ARTEFACTS
breakpoints <- breakpoints[type != "ARTEFACT"]
n_breaks = breakpoints[, .N]

# load hg38 ideogram
data(UCSC.HG38.Human.CytoBandIdeogram)

# manual aggregate by locus, since you can't melt them all together at once
breakpoints_melted = rbindlist(list(
    copy(breakpoints),
    copy(breakpoints)
))
#   copy relevant locus info from second position to first
breakpoints_melted[(n_breaks + 1):(2 * n_breaks), chr_from := chr_to]
breakpoints_melted[(n_breaks + 1):(2 * n_breaks), start_from := start_to]
breakpoints_melted[(n_breaks + 1):(2 * n_breaks), end_from := end_to]
breakpoints_melted[(n_breaks + 1):(2 * n_breaks), strand_from := strand_to]
#   drop second information since it's already duplicated
breakpoints_melted[, chr_to:= NULL]
breakpoints_melted[, start_to:= NULL]
breakpoints_melted[, end_to:= NULL]
breakpoints_melted[, strand_to:= NULL]
#   fix column names
breakpoints_melted <- breakpoints_melted[, .SD, .SDcols = c(
    "Sample", "chr_from", "start_from", "end_from", "strand_from", "score", "resolution"
)]
colnames(breakpoints_melted) = c(
    "Sample", "Chrom", "Start", "End",
    "Strand", "Odds", "Resolution"
)
# convert Chrom to ordered factor
breakpoints_melted[, Chrom := factor(Chrom, levels = CHRS, ordered = TRUE)]


# ==============================================================================
# Analysis
# ==============================================================================
# convert breakpoints into Mbp bins and counts the rearrangements
breakpoints_melted[, StartBin := floor(Start / 10^6)]
breakpoints_melted[, EndBin := floor(End / 10^6)]

# convert breakpoints in genomic regions into bins
breakpoints_binned = rbindlist(lapply(
    1:breakpoints_melted[, .N],
    function(i) {
        print(i)
        start_bin = breakpoints_melted[i, StartBin]
        end_bin = breakpoints_melted[i, EndBin]
        dt = data.table(
            Sample = breakpoints_melted[i, Sample],
            Chrom = breakpoints_melted[i, Chrom],
            Bin = start_bin:end_bin,
            Count = 1
        )
        return(dt)
    }
))
# count all breakpoints in each bin
breakpoints_summed = breakpoints_binned[, sum(Count), by = c("Sample", "Chrom", "Bin")]
colnames(breakpoints_summed) = c("Sample", "Chrom", "Bin", "Count")

fwrite(
    breakpoints_summed[order(Sample, Chrom, Bin)],
    "Breakpoints/Default/breakpoints.binned.tsv",
    sep = "\t",
    col.names = TRUE
)

# ==============================================================================
# Plots
# ==============================================================================
# set core component for circos plots
RCircos.Set.Core.Components(
    UCSC.HG38.Human.CytoBandIdeogram,
    NULL,
    tracks.inside = 0,
    tracks.outside = 0
)
# plot circos plots for each patient
for (s in SAMPLES) {
    cat(s, "\n")
    png(
        paste0("Plots/", s, ".circos.png"),
        width = 12,
        height = 12,
        units = "cm",
        res = 300
    )
    RCircos.Set.Plot.Area()
    RCircos.Chromosome.Ideogram.Plot()
    RCircos.Link.Plot(
        link.data = breakpoints[
            Sample == s,
            .(chr_from, start_from, end_from, chr_to, start_to, end_to)
        ],
        track.num = 1,
        by.chromosome = TRUE
    )
    # change colour of inter vs intra chromosomal events
    # params = RCircos.Get.Plot.Parameters()
    # params$PlotColor = breakpoints$Colour
    # RCircos.Reset.Plot.Parameters(params)
    dev.off()
}

# number of SVs detected per patient
gg = (
    ggplot(data = breakpoints[, .N, by = Sample])
    + geom_col(aes(x = Sample, y = N, fill = Sample))
    + labs(x = NULL, y = "Breakpoint pairs")
    + guides(fill = FALSE)
    + theme_minimal()
    + theme(
        axis.text.x = element_text(angle = 90, vjust = 0, hjust = 0.5)
    )
)
ggsave(
    "Plots/sv-counts.png",
    height = 12,
    width = 20,
    units = "cm"
)

gg = (
    ggplot(data = breakpoints_melted[, .N, by = c("Sample", "Chrom")])
    + geom_col(aes(x = Sample, y = N, fill = Sample))
    + labs(x = NULL, y = "Number of SVs")
    + guides(fill = guide_legend(title = "Patient"))
    + facet_grid(. ~ Chrom)
    + theme_minimal()
    + theme(
        axis.text.x = element_blank(),
        legend.position = "bottom"
    )
)
ggsave(
    "Plots/sv-counts-by-chrom.png",
    height = 12,
    width = 40,
    units = "cm"
)
ggsave(
    "Plots/sv-counts-by-chrom.pdf",
    height = 12,
    width = 40,
    units = "cm"
)

for (s in SAMPLES) {
    cat(s, "\n")
    gg = (
        ggplot(data = breakpoints_melted[Sample == s, .N, by = "Chrom"])
        + geom_col(aes(x = Chrom, y = N))
        + labs(x = NULL, y = "Number of SVs")
        + guides(fill = FALSE)
        + scale_x_discrete(drop=FALSE)
        + ylim(0, 31)
        + theme_minimal()
        + theme(
            axis.text.x = element_text(angle = 90, vjust = 0, hjust = 0.5),
            panel.grid.major.x = element_blank()
        )
    )
    ggsave(
        paste0("Plots/", s, ".sv-counts-by-chrom.png"),
        height = 12,
        width = 40,
        units = "cm"
    )
    ggsave(
        paste0("Plots/", s, ".sv-counts-by-chrom.pdf"),
        height = 12,
        width = 40,
        units = "cm"
    )
}

# number and location of SVs across all patients
gg = (
    ggplot()
    + geom_col(
        data = breakpoints_summed,
        mapping = aes(x = Bin, y = Count, fill = Sample),
        position = "stack"
    )
    # add geom_blank to ensure each facet is the correct size
    + geom_blank(data = hg38, mapping = aes(x = Bins, y = 0))
    + labs(x = "Position", y = "Count")
    + facet_grid(. ~ Chrom, scales = "free_x", space = "free_x", switch = "x")
    + theme_classic()
    + theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing = unit(0, "lines"),
        strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5)
    )
)
ggsave(
    "Plots/sv-loci.png",
    height = 12,
    width = 40,
    units = "cm"
)

gg = (
    ggplot(data = breakpoints_summed[, length(unique(Sample)), by = c("Chrom", "Bin")])
    + geom_bar(aes(x = V1))
    + labs(x = "Number of patients", y = "Number of Mbp bins with a breakpoint")
    + scale_x_discrete(
        limits = seq(1, 13),
        breaks = seq(1, 13, by = 2),
        labels = seq(1, 13, by = 2)
    )
    + theme_minimal()
)
ggsave(
    "Plots/sv-recurrence.png",
    height = 12,
    width = 20,
    units = "cm"
)
