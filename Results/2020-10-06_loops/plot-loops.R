# ==============================================================================
# Meta
# ==============================================================================
# Plot loops
# --------------------------------------
# Description: Plot statistical information about loop calls
# Author: James Hawley


# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("logging"))

loginfo("Loading packages")
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
suppressMessages(library("pheatmap"))
suppressMessages(library("Matrix"))

# ==============================================================================
# Functions
# ==============================================================================

# ==============================================================================
# Data
# ==============================================================================
loginfo("Loading data")
loops <- fread("Loops/merged-loops.sample-counts.tsv")

loop_jaccard <- fread("Loops/by-sample/loop-multi-jaccard.csv", sep=",", drop = 1)
colnames(loop_jaccard) <- gsub(
    "^.\\\\",
    "",
    gsub(
        ".flat-loops.bed$",
        "",
        colnames(loop_jaccard)
    )
)


# ==============================================================================
# Analysis
# ==============================================================================
loginfo("Performing calculations")
# calculate distances between loop anchors
loops[, Distance := pmin(abs(start_x - end_y), abs(end_x - start_y))]

distance_stats <- data.table(
    Statistic = c("Mean", "Median", "Mode"),
    Value = c(
        loops[, mean(Distance)],
        loops[, median(Distance)],
        loops[, as.numeric(names(which.max(table(Distance))))]
    )
)

# construct a symmetric jaccard matrix out of loop comparisons
jaccard_mtx <- as.matrix(loop_jaccard)
rownames(jaccard_mtx) <- colnames(jaccard_mtx)
jaccard_mtx[is.na(jaccard_mtx)] <- 0
jaccard_mtx <- jaccard_mtx + t(jaccard_mtx)
diag(jaccard_mtx) <- diag(jaccard_mtx) / 2

# ==============================================================================
# Save data
# ==============================================================================
loginfo("Saving data")

# ==============================================================================
# Plots
# ==============================================================================
loginfo("Plotting data")
gg <- (
    ggplot(data = loops)
    + geom_density(aes(x = Distance / 1000))
    + scale_x_continuous(
        name = "Distance between loop anchors (kbp)"
    )
    + geom_vline(
        data = distance_stats,
        mapping = aes(xintercept = Value / 1000, colour = Statistic),
        linetype = "dashed"
    )
    + geom_text(
        data = distance_stats,
        mapping = aes(
            x = Value / 1000,
            y = 0,
            colour = Statistic,
            label = paste(Statistic, "=", round(Value / 1000))
        ),
        angle = 90,
        hjust = 0,
        vjust = 1.1
    )
    + scale_colour_manual(
        limits = c("Mean", "Median", "Mode"),
        values = c("#1e90ff", "#66cdaa", "#b22222")
    )
    + guides(colour = FALSE)
    + theme_minimal()
)
ggsave("Plots/loop-calls.anchor-distance.png", gg, width = 12, height = 8, units = "cm")

pheatmap(
    mat=jaccard_mtx,
    filename="Plots/loop-jaccard.png"
)
