# ==============================================================================
# Meta
# ==============================================================================
# plot-grn-exprs
# --------------------------------------
# Description: Plot the expression changes for different types of GRNs
# Author: James Hawley


# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("logging"))

loginfo("Loading packages")
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))

# ==============================================================================
# Functions
# ==============================================================================
#' Save figures in multiple formats
#'
#' @param gg ggplot object
#' @param prefix Prefix for output file
#' @param ext Output extensions
#' @param dpi DPI resolution
savefig <- function(gg, prefix, ext = c("png", "pdf"), width = 20, height = 12, dpi = 400) {
    for (e in ext) {
        ggsave(
            paste(prefix, e, sep = "."),
            gg,
            height = height,
            width = width,
            units = "cm",
            dpi = dpi
        )
    }
}

# ==============================================================================
# Data
# ==============================================================================
loginfo("Loading data")
grn_sat <- fread("Graphs/grn-satisfiability.tsv", sep="\t")

# load expression data
exprs_genes <- fread(
    "../2020-06-18_sv-disruption-expression/results.genes.tsv",
    sep="\t"
)
exprs_tx <- fread(
    "../2020-06-18_sv-disruption-expression/results.transcripts.tsv",
    sep="\t"
)

grn_stats <- fread("Graphs/grn-stats.tsv", sep = "\t")

# ==============================================================================
# Analysis
# ==============================================================================
loginfo("Performing calculations")

# calculate summarize gene expression changes for each gene
exprs_summ <- exprs_tx[,
    .(
        n_transcripts = .N,
        Median_log2FC = median(b),
        Mean_log2FC = mean(b),
        WeightedMean_log2FC = sum(mean_obs * b) / sum(mean_obs)
    ),
    by = "ens_gene"
]

# merge GRN properties with gene expression changes
exprs <- merge(
    x = grn_sat,
    y = exprs_summ,
    by.x = "gene_id",
    by.y = "ens_gene",
    all.x = TRUE,
    all.y = FALSE
)

# ignore the GRNs that have no loops or no enhancers, since we can't predict anything from those
"%ni%" <- Negate("%in%")
exprs <- exprs[Rejection_Reason %ni% c("No loops", "No enhancers")]

fwrite(exprs, "Graphs/grn-exprs.tsv", sep = "\t")

# hypothesis testing
grns_only_gained_loops <- grn_stats[
    loops_gained > 0 & loops_shared >= 0 & loops_lost == 0 & enhancers_gained == 0 & enhancers_shared > 0 & enhancers_lost == 0,
    gene_id
]
grns_only_lost_loops <- grn_stats[
    loops_gained == 0 & loops_shared >= 0 & loops_lost > 0 & enhancers_gained == 0 & enhancers_shared > 0 & enhancers_lost == 0,
    gene_id
]
grns_only_gained_enhns <- grn_stats[
    loops_gained == 0 & loops_shared > 0 & loops_lost == 0 & enhancers_gained > 0 & enhancers_shared >= 0 & enhancers_lost == 0,
    gene_id
]
grns_only_lost_enhns <- grn_stats[
    loops_gained == 0 & loops_shared > 0 & loops_lost == 0 & enhancers_gained == 0 & enhancers_shared >= 0 & enhancers_lost > 0,
    gene_id
]
grns_no_changes <- grn_stats[
    loops_gained == 0 & loops_shared > 0 & loops_lost == 0 & enhancers_gained == 0 & enhancers_shared > 0 & enhancers_lost == 0,
    gene_id
]

t.test(
    x = grn_stats[gene_id %in% grns_only_gained_loops, Mean_log2FC],
    y = grn_stats[gene_id %in% grns_no_changes, Mean_log2FC],
    alternative = "greater"
)
t.test(
    x = grn_stats[gene_id %in% grns_only_gained_enhns, Mean_log2FC],
    y = grn_stats[gene_id %in% grns_no_changes, Mean_log2FC],
    alternative = "greater"
)
t.test(
    x = grn_stats[gene_id %in% grns_only_lost_loops, Mean_log2FC],
    y = grn_stats[gene_id %in% grns_no_changes, Mean_log2FC],
    alternative = "less"
)
t.test(
    x = grn_stats[gene_id %in% grns_only_lost_enhns, Mean_log2FC],
    y = grn_stats[gene_id %in% grns_no_changes, Mean_log2FC],
    alternative = "less"
)

# conert GRN stats to long format
grn_stats_long <- melt(
    grn_stats,
    id.vars = "gene_id",
    variable.name = "feature",
    value.name = "N"
)

# add fold change in expression to GRN stats
grn_stats <- merge(
    x = grn_stats,
    y = exprs[, .SD, .SDcols = c("gene_id", "Mean_log2FC")],
    by = "gene_id",
    all.x = TRUE
)

# grn_stats_model <- glm(
#     # Mean_log2FC ~ log10(loops_gained + 1) + log10(loops_shared + 1) + log10(loops_lost + 1) + log10(enhancers_gained + 1) + log10(enhancers_shared + 1) + log10(enhancers_lost + 1),
#     formula = Mean_log2FC ~ .,
#     dat = grn_stats,
#     family = "binomial"
# )


# ==============================================================================
# Plots
# ==============================================================================
loginfo("Plotting data")

exprs[, GRN_Class := gsub("^\\(", "", GRN_Class)]
exprs[, GRN_Class := gsub("\\)$", "", GRN_Class)]
exprs[, GRN_Class := gsub(", ", "\n", GRN_Class)]

text_y = 3.2

gg <- (
    ggplot(
        data = exprs,
        mapping = aes(
            x = GRN_Class,
            y = Mean_log2FC,
            group = GRN_Class
        )
    )
    + geom_boxplot(
        outlier.shape = NA, width=0.3, alpha = 0.1
    )
    + geom_point(
        aes(colour = GRN_Class),
        alpha = 0.2, position = position_jitter(height = 0, width = 0.3)
    )
    + geom_text(
        data = exprs[complete.cases(exprs), .N, by = GRN_Class],
        mapping = aes(
            x = GRN_Class,
            y = text_y,
            label = paste0("(", N, ")")
        ),
        vjust = -0.2
    )
    + scale_x_discrete(
        name = "GRN Alteration in T2E+"
    )
    + labs(y = bquote(log[2] * "(T2E+ / T2E-) expression"))
    + guides(fill = FALSE, colour = FALSE)
    + theme_minimal()
    + theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    )
)
savefig(gg, "Plots/grn-exprs.mean")

gg_grn_stats <- (
    ggplot(data = grn_stats_long)
    + geom_histogram(aes(x = N + 1), bins = sqrt(grn_stats[, .N]))
    + labs(x = "Count + 1", y = "Density")
    + facet_wrap(~ feature)
    + scale_x_log10()
    + theme_minimal()
)
savefig(gg_grn_stats, "Plots/grn-stats.dist")

# gg_grn_stats_fit <- (
#     ggplot()
#     + geom_point(
#         aes(
#             x = grn_stats[!is.na(Mean_log2FC), Mean_log2FC],
#             y = grn_stats_model$fitted.values
#         )
#     )
#     + labs(x = bquote("Observed " * log[2] * "(FC)"), y = bquote("Predicted " * log[2] * "(FC)"))
#     + theme_minimal()
# )
# savefig(gg_grn_stats_fit, "Plots/grn-stats.fit")

# gg_grn_stats_residuals <- (
#     ggplot()
#     + geom_point(
#         aes(
#             x = grn_stats[!is.na(Mean_log2FC), Mean_log2FC],
#             y = grn_stats_model$residuals
#         )
#     )
#     + labs(x = bquote("Observed " * log[2] * "(FC)"), y = "Residuals")
#     + theme_minimal()
# )
# savefig(gg_grn_stats_residuals, "Plots/grn-stats.fit.residuals")
