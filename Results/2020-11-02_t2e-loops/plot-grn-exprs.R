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

grn_stats <- merge(
    x = grn_stats,
    y = grn_sat[, .SD, .SDcols = c("gene_id", "GRN_Class")],
    by = "gene_id"
)

events <- list(
    "changed loop given changed enhancer" = grn_stats[
        enhancers_gained + enhancers_lost > 0,
        .(prob = .N / grn_stats[
            enhancers_gained + enhancers_lost > 0,
            .N
        ]),
        by = (loops_gained + loops_lost > 0)
    ],
    "changed enhancer given changed loop" = grn_stats[
        loops_gained + loops_lost > 0,
        .(prob = .N / grn_stats[
            loops_gained + loops_lost > 0,
            .N
        ]),
        by = (enhancers_gained + enhancers_lost > 0)
    ],
    "changed exprs given changed loop" = grn_stats[
        (loops_gained + loops_lost > 0) & !is.na(qval),
        .(prob = .N / grn_stats[
            (loops_gained + loops_lost > 0) & !is.na(qval),
            .N
        ]),
        by = (qval < 0.05)
    ],
    "changed exprs given changed enhancer" = grn_stats[
        (enhancers_gained + enhancers_lost > 0) & !is.na(qval),
        .(prob = .N / grn_stats[
            (enhancers_gained + enhancers_lost > 0) & !is.na(qval),
            .N
        ]),
        by = (qval < 0.05)
    ],
    "changed exprs given changed enhancer or changed loop" = grn_stats[
        ((loops_gained + loops_lost > 0) | (enhancers_gained + enhancers_lost > 0)) & !is.na(qval),
        .(prob = .N / grn_stats[
            ((loops_gained + loops_lost > 0) | (enhancers_gained + enhancers_lost > 0)) & !is.na(qval),
            .N
        ]),
        by = (qval < 0.05)
    ],
    "changed exprs given changed enhancer and changed loop" = grn_stats[
        (loops_gained + loops_lost > 0) & (enhancers_gained + enhancers_lost > 0) & !is.na(qval),
        .(prob = .N / grn_stats[
            (loops_gained + loops_lost > 0) & (enhancers_gained + enhancers_lost > 0) & !is.na(qval),
            .N
        ]),
        by = (qval < 0.05)
    ],
    "changed exprs given changed enhancer and same loop" = grn_stats[
        (loops_gained + loops_lost == 0) & (enhancers_gained + enhancers_lost > 0) & !is.na(qval),
        .(prob = .N / grn_stats[
            (loops_gained + loops_lost == 0) & (enhancers_gained + enhancers_lost > 0) & !is.na(qval),
            .N
        ]),
        by = (qval < 0.05)
    ],
    "changed exprs given same enhancer and changed loop" = grn_stats[
        (loops_gained + loops_lost > 0) & (enhancers_gained + enhancers_lost == 0) & !is.na(qval),
        .(prob = .N / grn_stats[
            (loops_gained + loops_lost > 0) & (enhancers_gained + enhancers_lost == 0) & !is.na(qval),
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
    "changed enhancer given changed exprs" = grn_stats[
        !is.na(qval) & (qval < 0.05),
        .(prob = .N / grn_stats[
            !is.na(qval) & (qval < 0.05),
            .N
        ]),
        by = (enhancers_gained + enhancers_lost > 0)
    ],
    "changed enhancer given changed exprs" = grn_stats[
        !is.na(qval) & (qval < 0.05),
        .(prob = .N / grn_stats[
            !is.na(qval) & (qval < 0.05),
            .N
        ]),
        by = (enhancers_gained + enhancers_lost > 0)
    ],
    "changed enhancer and changed loop given changed exprs" = grn_stats[
        !is.na(qval) & (qval < 0.05),
        .(prob = .N / grn_stats[
            !is.na(qval) & (qval < 0.05),
            .N
        ]),
        by = ((enhancers_gained + enhancers_lost > 0) & (loops_gained + loops_lost > 0))
    ],
    "changed enhancer and same loop given changed exprs" = grn_stats[
        !is.na(qval) & (qval < 0.05),
        .(prob = .N / grn_stats[
            !is.na(qval) & (qval < 0.05),
            .N
        ]),
        by = ((enhancers_gained + enhancers_lost > 0) & (loops_gained + loops_lost == 0))
    ],
    "same enhancer and changed loop given changed exprs" = grn_stats[
        !is.na(qval) & (qval < 0.05),
        .(prob = .N / grn_stats[
            !is.na(qval) & (qval < 0.05),
            .N
        ]),
        by = ((enhancers_gained + enhancers_lost == 0) & (loops_gained + loops_lost > 0))
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

fwrite(grn_stats, "Graphs/grn-exprs.tsv", sep = "\t")

# # hypothesis testing
# t.test(
#     x = exprs[GRN_Class == "(Same loop(s), Lost enhancer(s))", Mean_log2FC],
#     y = exprs[GRN_Class == "(Same loop(s), Same enhancer(s))", Mean_log2FC],
#     alternative = "less"
# )
# t.test(
#     x = exprs[GRN_Class == "(Same loop(s), Gained enhancer(s))", Mean_log2FC],
#     y = exprs[GRN_Class == "(Same loop(s), Same enhancer(s))", Mean_log2FC],
#     alternative = "greater"
# )


# ==============================================================================
# Plots
# ==============================================================================
loginfo("Plotting data")

grn_stats[, GRN_Class := gsub("^\\(", "", GRN_Class)]
grn_stats[, GRN_Class := gsub("\\)$", "", GRN_Class)]
grn_stats[, GRN_Class := gsub(", ", "\n", GRN_Class)]

text_y = 3.2

gg <- (
    ggplot(
        data = grn_stats,
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
        data = grn_stats[,
            .N,
            by = GRN_Class
        ],
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

gg <- (
    ggplot(
        data = grn_stats[!is.na(qval) & (qval < 0.05)],
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
        data = grn_stats[
            !is.na(qval) & (qval < 0.05),
            .N,
            by = GRN_Class
        ],
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
savefig(gg, "Plots/grn-exprs.sig.mean")
