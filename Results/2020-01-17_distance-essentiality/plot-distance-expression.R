# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))

# ==============================================================================
# Constants
# ==============================================================================
metadata = fread("../../Data/External/LowC_Samples_Data_Available.tsv")
SAMPLES = paste0("PCa", metadata[, get("Sample ID")])

# ==============================================================================
# Data
# ==============================================================================
# load TSS distances for each patient
distances = rbindlist(lapply(
    SAMPLES,
    function(s) {
        dt = fread(paste0("Closest/", s, ".distance-dependency.tsv"))
        dt[, SampleID := s]
        return(dt)
    }
))

# load gene expression information for each patient
expression = fread("../../Data/External/CPC-GENE/CPC-GENE_Chen-2019_RNAseq_rsem_gene_FPKM.13-Low-C-only.tsv")
# melt into long form
expression_long = melt(
    expression,
    id.vars = c("EnsemblID", "Symbol"),
    variable.name = "SampleID",
    value.name = "Expression"
)


# ==============================================================================
# Analysis
# ==============================================================================
# calculate the quantile each gene belongs to, wrt the individual
# because a large percentage of genes are not expressed (i.e. reads == 0 FPKM)
# sort those separately from expressed genes
expression_long[Expression == 0, Level := "0"]
for (s in SAMPLES) {
    breaks = expression_long[SampleID == s & Expression > 0, quantile(Expression, 1:5 / 5, na.rm = TRUE)]
    breaks = c(0, breaks)
    expression_long[
        SampleID == s & Expression > 0,
        Level := cut(
            Expression,
            breaks = breaks,
            labels = as.character(1:5)
        )
    ]
}

# merge distance and expression information
# remove everything after "." in Ensembl IDs, since these may conflict and not merge
expression_long[, shortID := gsub("\\..+$", "", EnsemblID)]
distances[, shortID := gsub("\\..+$", "", Ensembl_ID)]
# merges 19535 genes, discards 391
combined = merge(
    x = expression_long,
    y = distances,
    by = c("shortID", "SampleID")
)

# ==============================================================================
# Plots
# ==============================================================================
gg = (
    ggplot(data = combined)
    + geom_density(aes(x = Fraction, colour = Level))
    + labs(y = "Density")
    + scale_x_continuous(
        limits = c(0, 0.5),
        breaks = 0:5 / 10,
        name = "TSS distance from boundary",
        labels = paste0(0:5 * 10, "%")
    )
    + scale_colour_manual(
        breaks = 0:5,
        name = "Expression Quantile",
        labels = c(
            "Not Expressed",
            "[0%, 20%)",
            "[20%, 40%)",
            "[40%, 60%)",
            "[60%, 80%)",
            "[80%, 100%]"
        ),
        values = c(
            "#000000",
            "#edf8fb",
            "#b2e2e2",
            "#66c2a4",
            "#2ca25f",
            "#006d2c"
        )
    )
    + facet_wrap(~ SampleID)
    + theme_minimal()
    + theme(
        axis.text.x = element_text(angle = 90)
    )
)
ggsave(
    "Plots/distance-density-by-expression.png",
    height = 12,
    width = 20,
    units = "cm"
)

# same as before, but an empirical CDF plot instead of KDE
gg = (
    ggplot(data = combined[complete.cases(combined)])
    + stat_ecdf(aes(x = Fraction, colour = Level), geom = "step")
    # + geom_density(aes(x = Mean_Fraction, colour = Level))
    + labs(y = "Cumulative Density")
    + scale_x_continuous(
        limits = c(0, 0.5),
        breaks = 0:5 / 10,
        name = "TSS distance from boundary",
        labels = paste0(0:5 * 10, "%")
    )
    + scale_colour_manual(
        breaks = 0:5,
        name = "Expression Quantile",
        labels = c(
            "Not Expressed",
            "[0%, 20%)",
            "[20%, 40%)",
            "[40%, 60%)",
            "[60%, 80%)",
            "[80%, 100%]"
        ),
        values = c(
            "#000000",
            "#edf8fb",
            "#b2e2e2",
            "#66c2a4",
            "#2ca25f",
            "#006d2c"
        )
    )
    + facet_wrap(~ SampleID)
    + theme_minimal()
    + theme(
        axis.text.x = element_text(angle = 90)
    )
)
ggsave(
    "Plots/distance-density-by-expression.ecdf.png",
    height = 12,
    width = 20,
    units = "cm"
)

# plot size of TADs, stratified by their essentiality
gg = (
    ggplot(data = combined)
    + geom_density(aes(x = end_TAD - start_TAD, colour = Level))
    + labs(y = "Density")
    + scale_x_log10(name = "Parent TAD size")
    + scale_colour_manual(
        breaks = 0:5,
        name = "Expression Quantile",
        labels = c(
            "Not Expressed",
            "[0%, 20%)",
            "[20%, 40%)",
            "[40%, 60%)",
            "[60%, 80%)",
            "[80%, 100%]"
        ),
        values = c(
            "#000000",
            "#edf8fb",
            "#b2e2e2",
            "#66c2a4",
            "#2ca25f",
            "#006d2c"
        )
    )
    + facet_wrap(~ SampleID)
    + theme_minimal()
    + theme(
        axis.text.x = element_text(angle = 90)
    )
)
ggsave(
    "Plots/-by-expression.png",
    height = 12,
    width = 20,
    units = "cm"
)

# same as before, but an empirical CDF plot instead of KDE
gg = (
    ggplot(data = combined[complete.cases(combined)])
    + stat_ecdf(aes(x = end_TAD - start_TAD, colour = Level), geom = "step")
    + labs(y = "Cumulative Density")
    + scale_x_log10(name = "Parent TAD size")
    + scale_colour_manual(
        breaks = 0:5,
        name = "Expression Quantile",
        labels = c(
            "Not Expressed",
            "[0%, 20%)",
            "[20%, 40%)",
            "[40%, 60%)",
            "[60%, 80%)",
            "[80%, 100%]"
        ),
        values = c(
            "#000000",
            "#edf8fb",
            "#b2e2e2",
            "#66c2a4",
            "#2ca25f",
            "#006d2c"
        )
    )
    + facet_wrap(~ SampleID)
    + theme_minimal()
    + theme(
        axis.text.x = element_text(angle = 90)
    )
)
ggsave(
    "Plots/TAD-size-density-by-expression.ecdf.png",
    height = 12,
    width = 20,
    units = "cm"
)
