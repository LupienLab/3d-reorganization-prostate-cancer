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

callout_names = c("AR", "ETV1", "FOXA1", "HOXB13", "MYC", "SOX9")
callouts = genes[name %in% callout_names]
callouts_melted = genes_melted[name %in% callout_names]

# ==============================================================================
# Plots
# ==============================================================================
gg = (
    ggplot(data = genes[complete.cases(LNCaP)])
    + geom_point(aes(x = Fraction, y = LNCaP), alpha = 0.1)
    + geom_smooth(aes(x = Fraction, y = LNCaP), method = "loess")
    + labs(y = "Essentiality in LNCaP")
    + annotate(
        geom = "text",
        x = callouts[, Fraction],
        y = callouts[, LNCaP],
        label = callouts[, name]
    )
    + scale_x_continuous(
        limits = c(0, 0.5),
        breaks = 0:5 / 10,
        name = "TSS distance from boundary",
        labels = paste(0:5 * 10, "%")
    )
    + theme_minimal()
)
ggsave(
    paste0(ARGS$prefix, ".proximity.png"),
    height = 12,
    width = 20,
    units = "cm"
)

gg = (
    ggplot(data = genes[complete.cases(LNCaP)])
    + geom_density2d(aes(x = Fraction, y = LNCaP))
    + labs(y = "Essentiality in LNCaP")
    # + annotate(
    #     geom = "text",
    #     x = genes[name %in% c("FOXA1", "MYC", "AR"), Fraction],
    #     y = genes[name %in% c("FOXA1", "MYC", "AR"), LNCaP],
    #     label = c("FOXA1", "MYC", "AR")
    # )
    + scale_x_continuous(
        limits = c(0, 0.5),
        breaks = 0:5 / 10,
        name = "TSS distance from boundary",
        labels = paste(0:5 * 10, "%")
    )
    + theme_minimal()
)
ggsave(
    paste0(ARGS$prefix, ".proximity-density.png"),
    height = 12,
    width = 20,
    units = "cm"
)

# plot essentiality across the TAD
gg = (
    ggplot(data = genes_melted[complete.cases(genes_melted)])
    + geom_point(aes(x = Fraction, y = Essentiality), alpha = 0.1)
    + geom_smooth(aes(x = Fraction, y = Essentiality), method = "loess")
    + labs(x = "TSS distance from boundary", y = "Essentiality")
    + geom_text(
        data = callouts_melted,
        mapping = aes(
            x = Fraction,
            y = Essentiality,
            label = name
        ),
        inherit.aes = FALSE
    )
    + scale_x_continuous(
        limits = c(0, 0.5),
        breaks = 0:5 / 10,
        name = "TSS distance from boundary",
        labels = paste(0:5 * 10, "%")
    )
    + facet_wrap(~ Cell)
    + theme_minimal()
)
ggsave(
    paste0(ARGS$prefix, ".proximity.all-lines.png"),
    height = 12,
    width = 20,
    units = "cm"
)

# for a given quintile of essentiality, what distribution do they have across a TAD?
gg = (
    ggplot(data = genes_melted[complete.cases(genes_melted)])
    # + geom_density(aes(x = Fraction, fill = Quintile, colour = Quintile), alpha = 0)
    + geom_histogram(aes(x = Fraction), bins = 50)
    + labs(x = "TSS distance from boundary", y = "Essentiality")
    + scale_x_continuous(
        limits = c(0, 0.5),
        breaks = 0:5 / 10,
        name = "TSS distance from boundary",
        labels = paste(0:5 * 10, "%")
    )
    + facet_wrap(~ Cell)
    + theme_minimal()
)
ggsave(
    paste0(ARGS$prefix, ".proximity-by-essentiality-quintile.png"),
    height = 12,
    width = 20,
    units = "cm"
)

## qq plot of gene distribution throughout a TAD vs uniform distribution
gg = (
    ggplot(
        data = genes,
        mapping = aes(sample = Fraction)
    )
    + geom_qq(distribution = qunif, alpha = 0.1)
    + scale_x_continuous(
        limits = c(0, 1),
        name = "Uniform distribution quantile"
    )
    + scale_y_continuous(
        limits = c(0, 0.5),
        breaks = 0:5 / 10,
        name = "TSS distance from boundary quantile"
    )
    + theme_minimal()
)
ggsave(
    paste0(ARGS$prefix, ".proximity.qqplot.png"),
    height = 12,
    width = 20,
    units = "cm"
)
