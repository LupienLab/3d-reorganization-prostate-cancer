# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
suppressMessages(library("ggrepel"))
source(
    file.path("..", "2020-02-19_chromoplexy", "plotting-helper.R")
)

RES_DIR <- file.path("..", "..", "results", "2020-06-12_sv-disruption-acetylation")
PLOT_DIR <- file.path(RES_DIR, "Plots", "exprs-acetyl")

QVAL_THRESH <- 0.05
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
delta <- fread(
    file.path(RES_DIR, "exprs-acetyl.tsv"),
    sep = "\t",
    header = TRUE
)
delta[, col := ifelse(
    !sig,
    "NS",
    ifelse(
        abs(log2acetyl_fc) < 0.1,
        "NS",
        ifelse(log2exprs_fc > 0, "Up", "Down")
    )
)]

# ==============================================================================
# Plots
# ==============================================================================
h_cor <- cor.test(
    x = delta$log2exprs_fc,
    y = delta$log2acetyl_fc,
    method = "spearman"
)

# ==============================================================================
# Plots
# ==============================================================================
gg_cor <- (
    ggplot()
    + geom_vline(aes(xintercept = FC_THRESH), linetype = "dashed")
    + geom_vline(aes(xintercept = -FC_THRESH), linetype = "dashed")
    + geom_point(
        data = delta[sig == FALSE],
        mapping = aes(x = log2exprs_fc, y = log2acetyl_fc, fill = col),
        shape = 21,
        colour = "#000000",
        size = 2,
        alpha = 0.5
    )
    + geom_point(
        data = delta[sig == TRUE],
        mapping = aes(x = log2exprs_fc, y = log2acetyl_fc, fill = col),
        shape = 21,
        colour = "#000000",
        size = 2,
        alpha = 1
    )
    + geom_label_repel(
        data = delta[sig == TRUE & abs(log2acetyl_fc) > 0.1],
        mapping = aes(x = log2exprs_fc, y = log2acetyl_fc, label = gene_name)
    )
    + annotate(
        geom = "text",
        x = 4,
        y = -4,
        label = bquote(rho * " = " * .(round(h_cor$estimate, 4)))
    )
    + scale_fill_manual(
        breaks = c(
            "Down",
            "NS",
            "Up"
        ),
        labels = c(
            "Under-expressed",
            "N.S.",
            "Over-expressed"
        ),
        values = c(
            "#0000cd",
            "#b9b9b9",
            "#ff6347"
        ),
        name = expression(FDR < 0.05)
    )
    + labs(
        x = expression(log[2] * "(Expression mut fold change)"),
        y = expression(log[2] * "(H3K27ac mut fold change)")
    )
    + theme_minimal()
    + theme(legend.position = "bottom")
)
savefig(gg_cor, file.path(PLOT_DIR, "genes"))
