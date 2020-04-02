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
breakpoints <- fread("Graphs/sv-breakpoints.tsv", sep = "\t", header = TRUE)
# convert chr to ordered factor
breakpoints[, chr := factor(chr, levels = CHRS, ordered = TRUE)]
breakpoint_pairs <- fread("Graphs/sv-breakpoints.paired.tsv", sep = "\t", header = TRUE)

# load hg38 ideogram
data(UCSC.HG38.Human.CytoBandIdeogram)

# ==============================================================================
# Analysis
# ==============================================================================
# convert breakpoints into Mbp bins and counts the rearrangements
breakpoints[, StartBin := floor(start / 10^6)]
breakpoints[, EndBin := floor(end / 10^6)]

# convert breakpoints in genomic regions into bins
breakpoints_binned = rbindlist(lapply(
    1:breakpoints[, .N],
    function(i) {
        start_bin = breakpoints[i, StartBin]
        end_bin = breakpoints[i, EndBin]
        dt = data.table(
            SampleID = breakpoints[i, SampleID],
            chr = breakpoints[i, chr],
            Bin = start_bin:end_bin,
            Count = 1
        )
        return(dt)
    }
))
# count all breakpoints in each bin
breakpoints_summed = breakpoints_binned[, sum(Count), by = c("SampleID", "chr", "Bin")]
colnames(breakpoints_summed) = c("SampleID", "chr", "Bin", "Count")

fwrite(
    breakpoints_summed[order(SampleID, chr, Bin)],
    "Statistics/breakpoints.binned.tsv",
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
        paste0("Plots/Circos/", s, ".circos.png"),
        width = 12,
        height = 12,
        units = "cm",
        res = 400
    )
    RCircos.Set.Plot.Area()
    RCircos.Chromosome.Ideogram.Plot()
    RCircos.Link.Plot(
        link.data = breakpoint_pairs[
            SampleID == s,
            .(chr_x, start_x, end_x, chr_y, start_y, end_y)
        ],
        track.num = 1,
        by.chromosome = TRUE
    )
    # change colour of inter vs intra chromosomal events
    # params = RCircos.Get.Plot.Parameters()
    # params$PlotColor = breakpoints$Colour
    # RCircos.Reset.Plot.Parameters(params)
    dev.off()
    png(
        paste0("Plots/Circos/", s, ".circos.pdf"),
        width = 12,
        height = 12,
        units = "cm",
        res = 400
    )
    RCircos.Set.Plot.Area()
    RCircos.Chromosome.Ideogram.Plot()
    RCircos.Link.Plot(
        link.data = breakpoint_pairs[
            SampleID == s,
            .(chr_x, start_x, end_x, chr_y, start_y, end_y)
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
    ggplot(data = breakpoints[, .N, by = SampleID])
    + geom_col(aes(x = SampleID, y = N, fill = SampleID))
    + labs(x = NULL, y = "Unique Breakpoints")
    + guides(fill = FALSE)
    + theme_minimal()
    + theme(
        axis.text.x = element_text(angle = 90, vjust = 0, hjust = 0.5)
    )
)
ggsave(
    "Plots/breakpoint-stats/sv-counts.png",
    height = 12,
    width = 20,
    units = "cm",
    dpi = 400
)

gg = (
    ggplot(data = breakpoints[, .N, by = c("SampleID", "chr")])
    + geom_col(aes(x = SampleID, y = N, fill = SampleID))
    + labs(x = NULL, y = "Breakpoints")
    + guides(fill = guide_legend(title = "Patient"))
    + facet_grid(. ~ chr)
    + theme_minimal()
    + theme(
        panel.grid.major.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "bottom"
    )
)
ggsave(
    "Plots/breakpoint-stats/sv-counts-by-chrom.labelled.png",
    height = 12,
    width = 40,
    units = "cm",
    dpi = 400
)
ggsave(
    "Plots/breakpoint-stats/sv-counts-by-chrom.labelled.pdf",
    height = 12,
    width = 40,
    units = "cm",
    dpi = 400
)

# save as above without any text
gg = (
    ggplot(data = breakpoints[, .N, by = c("SampleID", "chr")])
    + geom_col(aes(x = SampleID, y = N, fill = SampleID))
    + labs(x = NULL, y = NULL)
    + guides(fill = guide_legend(title = NULL))
    + facet_grid(. ~ chr)
    + theme_minimal()
    + theme(
        panel.grid.major.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "bottom",
        strip.text.x = element_blank()
    )
)
ggsave(
    "Plots/breakpoint-stats/sv-counts-by-chrom.png",
    height = 12,
    width = 40,
    units = "cm",
    dpi = 400
)
ggsave(
    "Plots/breakpoint-stats/sv-counts-by-chrom.pdf",
    height = 12,
    width = 40,
    units = "cm",
    dpi = 400
)

# number and location of SVs across all patients
gg = (
    ggplot()
    + geom_col(
        data = breakpoints_summed,
        mapping = aes(x = Bin, y = Count, fill = SampleID),
        position = "stack"
    )
    # add geom_blank to ensure each facet is the correct size
    + geom_blank(data = hg38, mapping = aes(x = Bins, y = 0))
    + labs(x = "Position", y = "Count")
    + facet_grid(. ~ chr, scales = "free_x", space = "free_x", switch = "x")
    + theme_classic()
    + theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing = unit(0, "lines"),
        strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5)
    )
)
ggsave(
    "Plots/breakpoint-stats/sv-loci.png",
    height = 12,
    width = 40,
    units = "cm",
    dpi = 400
)
ggsave(
    "Plots/breakpoint-stats/sv-loci.pdf",
    height = 12,
    width = 40,
    units = "cm",
    dpi = 400
)

gg = (
    ggplot(data = breakpoints_summed[, length(unique(SampleID)), by = c("chr", "Bin")])
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
    "Plots/breakpoint-stats/sv-recurrence.png",
    height = 12,
    width = 20,
    units = "cm",
    dpi = 400
)
ggsave(
    "Plots/breakpoint-stats/sv-recurrence.pdf",
    height = 12,
    width = 20,
    units = "cm",
    dpi = 400
)