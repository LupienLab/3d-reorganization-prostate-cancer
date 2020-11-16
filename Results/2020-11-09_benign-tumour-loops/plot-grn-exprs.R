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

# 1. Probabilities of expression or GRN changes
# --------------------------------------
grn_stats <- merge(
    x = grn_stats,
    y = exprs_genes,
    by.x = "gene_id",
    by.y = "target_id",
    all.x = TRUE,
    all.y = FALSE
)

events <- list(
    "changed exprs" = grn_stats[
        !is.na(qval) & (qval < 0.05),
        .(prob = .N / grn_stats[
            !is.na(qval),
            .N
        ])
    ],
    "changed loop" = grn_stats[
        (loops_gained + loops_lost > 0),
        .(prob = .N / grn_stats[,
            .N
        ])
    ],
    "changed exprs given changed loop" = grn_stats[
        (loops_gained + loops_lost > 0) & !is.na(qval),
        .(prob = .N / grn_stats[
            (loops_gained + loops_lost > 0) & !is.na(qval),
            .N
        ]),
        by = (qval < 0.05)
    ],
    "changed exprs given same loop" = grn_stats[
        (loops_gained + loops_lost == 0) & !is.na(qval),
        .(prob = .N / grn_stats[
            (loops_gained + loops_lost == 0) & !is.na(qval),
            .N
        ]),
        by = (qval < 0.05)
    ],
    "changed loop given changed exprs" = grn_stats[
        !is.na(qval) & (qval < 0.05),
        .(prob = .N / grn_stats[
            !is.na(qval) & (qval < 0.05),
            .N
        ]),
        by = (loops_gained + loops_lost > 0)
    ],
    "changed loop given no changed exprs" = grn_stats[
        !is.na(qval) & (qval >= 0.05),
        .(prob = .N / grn_stats[
            !is.na(qval) & (qval >= 0.05),
            .N
        ]),
        by = (loops_gained + loops_lost > 0)
    ],
    "same loop given changed exprs" = grn_stats[
        !is.na(qval) & (qval < 0.05),
        .(prob = .N / grn_stats[
            !is.na(qval) & (qval < 0.05),
            .N
        ]),
        by = (loops_gained + loops_lost == 0)
    ],
    "same loop given no changed exprs" = grn_stats[
        !is.na(qval) & (qval >= 0.05),
        .(prob = .N / grn_stats[
            !is.na(qval) & (qval >= 0.05),
            .N
        ]),
        by = (loops_gained + loops_lost == 0)
    ]
)

# 2. Summarizing total gene expression
# --------------------------------------
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
grn_stats <- merge(
    x = grn_stats,
    y = exprs_summ[, .SD, .SDcols = c("ens_gene", "Mean_log2FC")],
    by.x = "gene_id",
    by.y = "ens_gene",
    all.x = TRUE,
    all.y = FALSE
)

grn_stats[,
    Loop_Changes := ifelse(
        loops_gained > 0,
        ifelse(loops_lost > 0, "Opposing loop changes", "Gained loop(s)"),
        ifelse(loops_lost > 0, "Lost loop(s)", "Same loop(s)")
    )
]

fwrite(grn_stats, "Graphs/grn-exprs.tsv", sep = "\t")

# 3. Hypothesis testing for gene expression changes and GRN changes
# --------------------------------------
t.test(
    x = grn_stats[Loop_Changes == "Gained loop(s)", Mean_log2FC],
    y = grn_stats[Loop_Changes == "Same loop(s)", Mean_log2FC],
    alternative = "greater"
)
t.test(
    x = grn_stats[Loop_Changes == "Lost loop(s)", Mean_log2FC],
    y = grn_stats[Loop_Changes == "Same loop(s)", Mean_log2FC],
    alternative = "less"
)

grn_associations <- data.table(
    "Changed_Exprs" = rep(c(TRUE, FALSE), each = 2),
    "Changed_GRN" = factor(
        rep(c("Loop(s)", "No change"), 2),
        ordered = TRUE,
        levels = c("No change", "Loop(s)")
    ),
    "Proportion" = c(
        events[["changed loop given changed exprs"]][loops_gained == TRUE, prob],
        events[["same loop given changed exprs"]][loops_gained == TRUE, prob],
        events[["changed loop given no changed exprs"]][loops_gained == TRUE, prob],
        events[["same loop given no changed exprs"]][loops_gained == TRUE, prob]
    )
)


# ==============================================================================
# Plots
# ==============================================================================
loginfo("Plotting data")

text_y = 3.2

gg <- (
    ggplot(
        data = grn_stats,
        mapping = aes(
            x = Loop_Changes,
            y = Mean_log2FC,
            group = Loop_Changes
        )
    )
    + geom_violin(
        aes(fill = Loop_Changes)
    )
    + geom_boxplot(
        outlier.shape = NA,
        width = 0.3,
        alpha = 0.1
    )
    + geom_text(
        data = grn_stats[,
            .N,
            by = Loop_Changes
        ],
        mapping = aes(
            x = Loop_Changes,
            y = text_y,
            label = paste0("(", N, ")")
        ),
        vjust = -0.2
    )
    + scale_x_discrete(
        name = "GRN Alteration in Tumour"
    )
    + labs(y = bquote(log[2] * "(Tumour / Benign) expression"))
    + guides(fill = FALSE, colour = FALSE)
    + theme_minimal()
    + theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    )
)
savefig(gg, "Plots/grn-exprs.mean")

gg <- (
    ggplot(
        data = grn_stats,
        mapping = aes(
            x = 1,
            y = Mean_log2FC,
            group = Loop_Changes
        )
    )
    + geom_violin(
        aes(fill = Loop_Changes)
    )
    + geom_boxplot(
        outlier.shape = NA,
        width = 0.3,
        alpha = 0.1
    )
        + geom_text(
        data = grn_stats[,
            .N,
            by = Loop_Changes
        ],
        mapping = aes(
            x = Loop_Changes,
            y = text_y,
            label = paste0("(", N, ")")
        ),
        vjust = -0.2
    )
    + scale_x_discrete(
        name = "GRN Alteration in Tumour"
    )
    + labs(y = bquote(log[2] * "(Tumour / Benign) expression"))
    + guides(fill = FALSE, colour = FALSE)
    + facet_grid(~ Loop_Changes, scales = "free_x")
    + theme_minimal()
    + theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    )
)
savefig(gg, "Plots/grn-exprs.mean.facetted")

gg <- (
    ggplot(
        data = grn_stats[!is.na(qval) & (qval < 0.05)],
        mapping = aes(
            x = Loop_Changes,
            y = Mean_log2FC,
            group = Loop_Changes
        )
    )
    + geom_violin(
        aes(fill = Loop_Changes)
    )
    + geom_boxplot(
        outlier.shape = NA, width = 0.3, alpha = 0.1
    )
    + geom_text(
        data = grn_stats[
            !is.na(qval) & (qval < 0.05),
            .N,
            by = Loop_Changes
        ],
        mapping = aes(
            x = Loop_Changes,
            y = text_y,
            label = paste0("(", N, ")")
        ),
        vjust = -0.2
    )
    + scale_x_discrete(
        name = "GRN Alteration in Tumour"
    )
    + labs(y = bquote(log[2] * "(Tumour / Benign) expression"))
    + guides(fill = FALSE, colour = FALSE)
    + theme_minimal()
    + theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    )
)
savefig(gg, "Plots/grn-exprs.sig.mean")

gg <- (
    ggplot(
        data = grn_stats[!is.na(qval) & (qval < 0.05)],
        mapping = aes(
            x = 1,
            y = Mean_log2FC,
            group = Loop_Changes
        )
    )
    + geom_violin(
        aes(fill = Loop_Changes)
    )
    + geom_boxplot(
        outlier.shape = NA, width = 0.3, alpha = 0.1
    )
    + geom_point(
        # aes(colour = GRN_Class),
        alpha = 0.2, position = position_jitter(height = 0, width = 0.3)
    )
    # + geom_text(
    #     data = grn_stats[,
    #         .N,
    #         by = GRN_Class
    #     ],
    #     mapping = aes(
    #         x = GRN_Class,
    #         y = text_y,
    #         label = paste0("(", N, ")")
    #     ),
    #     vjust = -0.2
    # )
    + scale_x_discrete(
        name = "GRN Alteration in Tumour"
    )
    + labs(y = bquote(log[2] * "(Tumour / Benign) expression"))
    + guides(fill = FALSE, colour = FALSE)
    + facet_grid(~ Loop_Changes)
    + theme_minimal()
    + theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    )
)
savefig(gg, "Plots/grn-exprs.sig.mean.facetted")

# conert GRN stats to long format
grn_stats_long <- melt(
    grn_stats,
    id.vars = "gene_id",
    measure.vars = c("loops_gained", "loops_shared", "loops_lost"),
    variable.name = "feature",
    value.name = "N"
)

gg_grn_stats <- (
    ggplot(data = grn_stats_long)
    + geom_histogram(aes(x = N), bins = sqrt(grn_stats[, .N]))
    + labs(x = "Count + 1", y = "Density")
    + facet_wrap(~ feature)
    # + scale_x_log10()
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
gg <- (
    ggplot(data = grn_associations)
    + geom_col(
        aes(x = Changed_Exprs, y = Proportion, fill = Changed_GRN),
        position = position_stack()
    )
    + labs(x = "Differential gene expression", y = "Proportion")
    + scale_fill_manual(
        name = "GRN Alteration in Tumour",
        breaks = c(
            "Loop(s)",
            "Loop(s) and Enhancer(s)",
            "Enhancer(s)",
            "No change"
        ),
        values = c(
            "#fa8072",
            "#c67ca2",
            "#1e90ff",
            "#c0c0c0"
        )
    )
    + theme_minimal()
)
savefig(gg, "Plots/grn-exprs.change-prop")
