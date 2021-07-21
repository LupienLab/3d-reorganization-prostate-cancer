# ==============================================================================
# Meta
# ==============================================================================
# plot-acetyl-loops
# --------------------------------------
# Description: Plot results from combining T2E+/- differential acetylation and loop calls
# Author: James Hawley


# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("logging"))

loginfo("Loading packages")
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
source(file.path("..", "src", "savefig.R"))

LOOP_TYPES = c("T2E-specific", "nonT2E-specific", "shared")
LOOP_LABELS = c("T2E+", "T2E-", "Shared")


# ==============================================================================
# Data
# ==============================================================================
loginfo("Loading data")
hits <- fread("loops.intersected-acetyl.tsv")
hits[, loop_Loop_Type := factor(
    loop_Loop_Type,
    levels = LOOP_TYPES
)]

hits_sig <- fread("loops.intersected-sig-acetyl.tsv")
hits_sig[, loop_Loop_Type := factor(
    loop_Loop_Type,
    levels = LOOP_TYPES
)]

# ==============================================================================
# Analysis
# ==============================================================================
loginfo("Performing calculations")

# ==============================================================================
# Plots
# ==============================================================================
loginfo("Generating plots")
gg <- (
    ggplot(data = hits)
    + geom_point(
        aes(x = loop_Loop_Type, y = peak_Fold, colour = loop_Loop_Type),
        alpha = 0.1,
        position = position_jitter(width = 0.2)
    )
    + geom_boxplot(
        aes(x = loop_Loop_Type, y = peak_Fold, colour = loop_Loop_Type),
        outlier.shape = NA,
        width = 0.5,
        alpha = 0.8
    )
    + labs(y = "log2(Acetylation Fold Change) T2E - Non-T2E")
    + scale_x_discrete(
        name = "Loop Type",
        breaks = LOOP_TYPES,
        labels = LOOP_LABELS
    )
    + theme_minimal()
)
savefig(gg, "Plots/diff-loop-diff-acetyl.all")

gg <- (
    ggplot(data = hits_sig)
    + geom_point(
        aes(x = loop_Loop_Type, y = peak_Fold, colour = loop_Loop_Type),
        alpha = 0.1,
        position = position_jitter(width = 0.2)
    )
    + geom_boxplot(
        aes(x = loop_Loop_Type, y = peak_Fold, colour = loop_Loop_Type),
        outlier.shape = NA,
        width = 0.5,
        alpha = 0.8
    )
    + labs(y = "log2(Acetylation Fold Change) T2E - Non-T2E")
    + scale_x_discrete(
        name = "Loop Type",
        breaks = LOOP_TYPES,
        labels = LOOP_LABELS
    )
    + theme_minimal()
)
savefig(gg, "Plots/diff-loop-diff-acetyl.sig")
