# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("pheatmap"))
suppressMessages(library("RColorBrewer"))
suppressMessages(library("ggplot2"))
suppressMessages(library("ggrepel"))
source(
    file.path("..", "2020-02-19_chromoplexy", "plotting-helper.R")
)

RES_DIR <- file.path("..", "..", "results", "2020-06-12_sv-disruption-acetylation")
ACETYL_DIR <- file.path(RES_DIR, "Acetylation")
TEST_DIR <- file.path(ACETYL_DIR, "Tests")
PLOT_DIR <- file.path(RES_DIR, "Plots")

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
metadata <- metadata[Include == "Yes"]
SAMPLES <- metadata[, SampleID]

test_IDs <- sort(sapply(
    list.files(TEST_DIR, "local.tsv"),
    function(fn) {
        as.numeric(gsub("\\.local\\.tsv", "", gsub("test_", "", fn)))
    }
))

acetyl <- rbindlist(lapply(
    test_IDs,
    function(tid) {
        dt <- fread(
            file.path(
                TEST_DIR,
                paste0("test_", tid, ".local.tsv")
            ),
            sep = "\t",
            header = TRUE
        )
        dt[, test_ID := tid]
        return(dt)
    }
))
acetyl[, col := paste(
    (FDR < QVAL_THRESH),
    ifelse(
        abs(Fold) > FC_THRESH,
        ifelse(Fold > FC_THRESH, "Up", "Down"),
        "FALSE"
    ),
    sep = "_"
)]


# ==============================================================================
# Plots
# ==============================================================================
# QQ plot of local p-values
gg_qq <- (
    ggplot(data = acetyl[order(p.value)])
    + geom_point(aes(x = -log10(ppoints(acetyl[, .N])), y = -log10(p.value)))
    + geom_abline(aes(intercept = 0, slope = 1), linetype = "dashed")
    + scale_x_continuous(
        limits = c(0, 3)
    )
    + scale_y_continuous(
        limits = c(0, 15)
    )
    + labs(
        x = expression(-log[10] * "(Uniform quantile)"),
        y = expression(-log[10] * "(Observed quantile)")
    )
    + theme_minimal()
)
savefig(gg_qq, file.path(PLOT_DIR, "Distribution", "acetylation.qq-plot"))

# histogram of p-values
gg_pval_hist <- (
    ggplot(data = acetyl)
    + geom_histogram(aes(x = p.value))
    + labs(x = "p-value", y = "Frequency")
    + theme_minimal()
)
savefig(gg_pval_hist, file.path(PLOT_DIR, "Distribution", "acetylation.p-values"))

# volcano plot of p-value vs log2 fold change
gg_volcano <- (
    ggplot(data = acetyl)
    + geom_vline(aes(xintercept = -LOG2FOLD_THRESH), linetype = "dashed")
    + geom_vline(aes(xintercept = LOG2FOLD_THRESH), linetype = "dashed")
    + geom_hline(aes(yintercept = -log10(QVAL_THRESH)), linetype = "dashed")
    + geom_point(aes(x = Fold, y = -log10(FDR), colour = col))
    + scale_colour_manual(
        breaks = c(
            "FALSE_Down",
            "FALSE_FALSE",
            "FALSE_Up",
            "TRUE_Down",
            "TRUE_FALSE",
            "TRUE_Up"
        ),
        values = c(
            "#777777",
            "#b9b9b9",
            "#777777",
            "#0000cd",
            "#777777",
            "#ff6347"
        )
    )
    + labs(x = expression(log[2] * "(Fold Change)"), y = expression(-log[10] * " FDR"))
    + guides(colour = FALSE)
    + theme_minimal()
)
savefig(gg_volcano, file.path(PLOT_DIR, "Distribution", "acetylation.volcano"), ext="png")
