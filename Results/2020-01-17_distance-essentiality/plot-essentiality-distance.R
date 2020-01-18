# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
suppressMessages(library("argparse"))

if (!interactive()) {
    PARSER <- argparse::ArgumentParser(
        description = "Plot the distribution of essentiality scores of genes across the length of TADs"
    )
    PARSER$add_argument(
        "data",
        type = "character",
        help = "Input data file"
    )
    PARSER$add_argument(
        "-p", "--prefix",
        type = "character",
        help = "Prefix for output image files"
    )
    ARGS <- PARSER$parse_args()
} else {
    ARGS = list(
        data = "essential-distance/PCa3023.distance-dependency.tsv",
        prefix = "output"
    )
}

# ==============================================================================
# Data
# ==============================================================================
cells = c("22Rv1", "DU145", "LNCaP", "MDA PCa 2b", "NCI-H660", "PC3", "VCaP")
# read in data
genes = fread(ARGS$data, sep = "\t", header = TRUE)

genes_melted = melt(
    genes,
    measure.vars = cells,
    variable.name = "Cell",
    value.name = "Essentiality"
)

for (cell in cells) {
    idx = genes_melted[, which(Cell == cell & !is.na(Essentiality))]
    breaks = genes_melted[idx, quantile(Essentiality, 0:5/5)]
    genes_melted[idx, Quintile := cut(Essentiality, breaks = breaks, labels = 1:5)]
}

# ==============================================================================
# Plots
# ==============================================================================
gg = (
    ggplot(data = genes[complete.cases(LNCaP)])
    + geom_smooth(aes(x = Fraction, y = LNCaP), method = "loess")
    + labs(x = "TSS distance from boundary", y = "Essentiality in LNCaP")
    + xlim(c(0, 0.5))
    + theme_minimal()
)
ggsave(
    paste0(ARGS$prefix, ".proximity.png"),
    height = 12,
    width = 20,
    units = "cm"
)

# plot essentiality across the TAD
gg = (
    ggplot(data = genes_melted[complete.cases(genes_melted)])
    + geom_smooth(aes(x = Fraction, y = Essentiality), method = "loess")
    + labs(x = "TSS distance from boundary", y = "Essentiality")
    + facet_wrap(~ Cell)
    + theme_minimal()
)
ggsave(
    paste0(ARGS$prefix, ".proximity.png"),
    height = 12,
    width = 20,
    units = "cm"
)

# for a given quintile of essentiality, what distribution do they have across a TAD?
gg = (
    ggplot(data = genes_melted[complete.cases(genes_melted)])
    + geom_density(aes(x = Fraction, fill = Quintile, colour = Quintile), alpha = 0)
    + labs(x = "TSS distance from boundary", y = "Essentiality")
    + xlim(c(0, 0.5))
    + facet_wrap(~ Cell)
    + theme_minimal()
)
ggsave(
    paste0(ARGS$prefix, ".proximity-by-essentiality-quintile.png"),
    height = 12,
    width = 20,
    units = "cm"
)
