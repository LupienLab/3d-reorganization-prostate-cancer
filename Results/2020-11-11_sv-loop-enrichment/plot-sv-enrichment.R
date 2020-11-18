# ==============================================================================
# Meta
# ==============================================================================
# plot-sv-enrichment
# --------------------------------------
# Description: Plot the distribution of permutation fold changes for overlaps of SVs that affect expression with loops
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
obs_enrichment <- fread("sv-loop-intersection.observed-enrichment.tsv", sep = "\t")
obs_log2fc <- (
    obs_enrichment[DGE == TRUE, log2(Median_Loops)]
    - obs_enrichment[DGE == FALSE, log2(Median_Loops)]
)
perm_fcs <- fread("sv-loop-intersection.permutations.tsv")

pval <- perm_fcs[dge_log2Fold >= obs_log2fc, .N] / perm_fcs[, .N]

loop_int <- fread("sv-loop-intersection.tsv", sep = "\t")

# ==============================================================================
# Analysis
# ==============================================================================
# calculate power of the permutation test
power_delta <- function(delta) {
    critical_val <- perm_fcs[, quantile(dge_log2Fold, 0.95)]
    return(1 - perm_fcs[dge_log2Fold <= critical_val - delta, .N] / perm_fcs[, .N])
}
cat("Power of |log2FC| = 1.0:", power_delta(1), "\n")
cat("Power of |log2FC| = 1.1:", power_delta(1.1), "\n")
cat("Power of |log2FC| = 1.2:", power_delta(1.2), "\n")
cat("Power of |log2FC| = 1.3:", power_delta(1.3), "\n")
cat("Power of |log2FC| = 1.4:", power_delta(1.4), "\n")
cat("Power of |log2FC| = 1.5:", power_delta(1.5), "\n")

# ==============================================================================
# Plots
# ==============================================================================
loginfo("Plotting permutations")
gg <- (
    ggplot(data = perm_fcs)
    + geom_density(
        aes(x = dge_log2Fold),
        fill = "#ba55d3",
        alpha = 0.2
    )
    + geom_vline(
        aes(xintercept = obs_log2fc),
        linetype = "dashed",
        colour = "#ba55d3"
    )
    + geom_text(
        aes(x = obs_log2fc, y = 0.6, label = "Observed fold change"),
        colour = "#ba55d3",
        angle = 90,
        hjust = 0.5,
        vjust = 1.3
    )
    + geom_text(
        aes(x = 1.5, y = 0.6, label = paste0("p = ", pval))
    )
    + labs(
        x = bquote(log[2] * "(Median loops FC)\n\n(Functional / Non-functional SVs)"),
        y = "Density"
    )
    + theme_minimal()
)
# subset region and plot
d <- ggplot_build(gg)$data[[1]]
gg <- gg + geom_area(
    data = subset(d, x >= obs_log2fc),
    aes(x = x, y = y),
    fill = "#ba55d3"
)
savefig(gg, "Plots/permutations")


gg_obs <- (
    ggplot(
        data = loop_int,
        mapping = aes(x = any_dge, y = n_loops)
    )
    + geom_point(
        aes(fill = any_dge),
        position = position_jitter(height = 0, width = 0.3),
        alpha = 0.7,
        shape = 21,
        size = 3
    )
    + geom_boxplot(colour = "#000000", alpha = 0.7, width = 0.2, outlier.shape = NA)
    + geom_path(
        data = data.table(
            x = c(FALSE, FALSE, TRUE, TRUE),
            y = c(125, 130, 130, 110),
            group = 1
        ),
        mapping = aes(x = x, y = y, group = group),
        linetype = "solid",
        colour = "#000000"
    )
    + annotate(
        geom = "text",
        x = 1.5,
        y = 130,
        label = paste0("p = ", pval),
        hjust = 0.5,
        vjust = -1,
        size = 4
    )
    + scale_x_discrete(
        name = "SV-Related Differential\nExpression",
        breaks = c(FALSE, TRUE),
        labels = c("No", "Yes")
    )
    + scale_y_continuous(
        name = "Loops",
        limits = c(0, 150)
    )
    + scale_fill_manual(
        breaks = c(FALSE, TRUE),
        values = c("#AEA28E", "#FC6347")
    )
    + guides(colour = FALSE, fill = FALSE)
    + theme_minimal()
)
savefig(gg_obs, "Plots/obs-loops", height = 6.5, width = 7)
