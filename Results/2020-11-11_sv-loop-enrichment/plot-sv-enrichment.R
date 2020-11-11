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

# ==============================================================================
# Data
# ==============================================================================
loginfo("Loading data")
obs <- fread("sv-loop-intersection.observed-enrichment.tsv", sep = "\t")
obs_fc <- obs[DGE == TRUE, Median_Loops] / obs[DGE == FALSE, Median_Loops]
perm_fcs <- fread("sv-loop-intersection.permutations.tsv")

pval <- perm_fcs[DGE_Fold >= obs_fc, .N] / perm_fcs[, .N]

# ==============================================================================
# Plots
# ==============================================================================
loginfo("Plotting permutations")
gg <- (
    ggplot(data = perm_fcs)
    + geom_density(
        aes(x = DGE_Fold),
        fill = "#ba55d3",
        alpha = 0.2
    )
    + geom_vline(
        aes(xintercept = obs_fc),
        linetype = "dashed",
        colour = "#ba55d3"
    )
    + geom_text(
        aes(x = obs_fc, y = 0.6, label = "Observed fold change"),
        colour = "#ba55d3",
        angle = 90,
        hjust = 0.5,
        vjust = 1.3
    )
    + geom_text(
        aes(x = 2.5, y = 0.6, label = paste0("p = ", pval))
    )
    + labs(x = "Median loop fold change (DGE / non-DGE)", y = "Density")
    + theme_minimal()
)
# subset region and plot
d <- ggplot_build(gg)$data[[1]]
gg <- gg + geom_area(
    data = subset(d, x >= obs_fc),
    aes(x = x, y = y),
    fill = "#ba55d3"
)

ggsave(
    "Plots/permutations.png",
    width = 20,
    height = 12,
    units = "cm"
)
