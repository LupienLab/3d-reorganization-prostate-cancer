# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
suppressMessages(library("ggrepel"))
source("../2020-02-19_chromoplexy/plotting-helper.R")

PLOT_DIR <- file.path("Plots", "exprs-acetyl")
FC_THRESH <- 1

# ==============================================================================
# Data
# ==============================================================================
# load sample metadata
metadata <- fread(
    "config.tsv",
    sep = "\t",
    header = TRUE
)
low_qual_samples <- metadata[Include == "No"]
metadata <- metadata[Include == "Yes"]
SAMPLES <- metadata$SampleID

# load comparison data
delta <- fread("exprs-acetyl.tsv", sep = "\t", header = TRUE)
delta[, col := paste((abs(acetyl) > FC_THRESH), (abs(exprs * log2(exp(1))) > FC_THRESH), sep = "_")]

# ==============================================================================
# Plots
# ==============================================================================
h_cor <- cor.test(
    x = delta$exprs,
    y = delta$acetyl,
    method = "spearman"
)

# ==============================================================================
# Plots
# ==============================================================================
gg_cor <- (
    ggplot()
    + geom_hline(aes(yintercept = FC_THRESH), linetype = "dashed")
    + geom_hline(aes(yintercept = -FC_THRESH), linetype = "dashed")
    + geom_vline(aes(xintercept = FC_THRESH), linetype = "dashed")
    + geom_vline(aes(xintercept = -FC_THRESH), linetype = "dashed")
    + geom_point(
        data = delta,
        mapping = aes(x = exprs * log2(exp(1)), y = acetyl, colour = col)
    )
    + geom_label_repel(
        data = delta[col == "TRUE_TRUE"],
        mapping = aes(x = exprs * log2(exp(1)), y = acetyl, label = gene_name)
    )
    + annotate(
        geom = "text",
        x = 4,
        y = -4,
        label = bquote(rho * " = " * .(round(h_cor$estimate, 4)))
    )
    + scale_colour_manual(
        breaks = c(
            "FALSE_FALSE",
            "FALSE_TRUE",
            "TRUE_FALSE",
            "TRUE_TRUE"
        ),
        values = c(
            "#b9b9b9",
            "#777777",
            "#3b3b3b",
            "#ff6347"
        )
    )
    + labs(x = expression(log[2] * "(Expression fold change)"), y = expression(log[2] * "(H3K27ac fold change)"))
    + guides(colour = FALSE)
    + xlim(-6, 6)
    + ylim(-6, 6)
    + theme_minimal()
)
savefig(gg_cor, file.path(PLOT_DIR, "genes"))

