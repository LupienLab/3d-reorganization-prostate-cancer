# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
suppressMessages(library("argparse"))
suppressMessages(library("RCircos"))

if (!interactive()) {
    PARSER <- argparse::ArgumentParser(
        description = "Make circos and other plots for breakpoint calls"
    )
    PARSER$add_argument(
        "breaks",
        type = "character",
        help = "Aggregated breakpoint call TSV file"
    )
    ARGS <- PARSER$parse_args()
} else {
    ARGS = list(
        breaks = "Breakpoints/Default/breakpoints.tsv"
    )
}

# ==============================================================================
# Functions
# ==============================================================================

# ==============================================================================
# Data
# ==============================================================================
# load breakpoint data
breakpoints = fread(ARGS$breaks, sep = "\t", header = TRUE)

# load hg38 ideogram
data(UCSC.HG38.Human.CytoBandIdeogram)

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
for (s in unique(breakpoints$Sample)) {
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
            .(Chrom1, Start1, End1, Chrom2, Start2, End2)
        ],
        track.num = 1,
        by.chromosome = FALSE
    )
    dev.off()
}

# number of SVs detected per patient
gg = (
    ggplot(data = breakpoints[, .N, by = Sample])
    + geom_col(aes(x = Sample, y = N, fill = Sample))
    + labs(x = "Patient", y = "Number of SVs")
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
