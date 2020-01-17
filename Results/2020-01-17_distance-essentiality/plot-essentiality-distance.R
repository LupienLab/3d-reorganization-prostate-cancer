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
        prefix = "Plots/PCa3023"
    )
}

# ==============================================================================
# Data
# ==============================================================================
# read in data
genes = fread(ARGS$data, sep = "\t", header = TRUE)
genes_melted = melt(
    genes,
    measure.vars = c("DU145", "LNCaP", "MDA PCa 2b", "NCI-H660", "PC3", "VCaP"),
    variable.name = "Cell",
    value.name = "Essentiality"
)

# ==============================================================================
# Plots
# ==============================================================================
gg = (
    ggplot(data = genes_melted[complete.cases(genes_melted)])
    + geom_smooth(aes(x = Fraction, y = Essentiality), method = "loess")
    + labs(x = "Fraction along TAD", y = "Essentiality")
    + facet_wrap(~ Cell)
    + theme_minimal()
)
ggsave(
    paste0(ARGS$prefix, ".proximity.png"),
    height = 12,
    width = 20,
    units = "cm"
)
