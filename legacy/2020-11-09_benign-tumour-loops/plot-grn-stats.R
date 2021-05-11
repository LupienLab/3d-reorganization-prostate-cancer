# ==============================================================================
# Meta
# ==============================================================================
# plot-grn-stats
# --------------------------------------
# Description: Make plots about the statistics of GRNs and their structure
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
grn_stats <- fread("Graphs/grn-stats.tsv", sep = "\t", header = TRUE)

grn_stats[, total_loops := loops_gained + loops_lost + loops_shared]

grn_stats_long <- melt(
    grn_stats,
    id.vars = "gene_id",
    variable.name = "feature",
    value.name = "N"
)

# ==============================================================================
# Plots
# ==============================================================================
loginfo("Making plots")
gg <- (
    ggplot(data = grn_stats)
    + geom_histogram(
        aes(x = total_loops),
        binwidth = 1
    )
    + labs(x = "Detected loops in a GRN", y = "Count")
    + theme_minimal()
)
savefig(gg, "Plots/grn-loops")
