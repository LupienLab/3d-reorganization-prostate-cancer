# ==============================================================================
# Meta
# ==============================================================================
# plot-shared-enhancers
# --------------------------------------
# Description: Plot the results of the shared enhancer permutation test
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
shared_enhns <- fread("shared-enhancers.tsv", sep = "\t")

# ==============================================================================
# Plots
# ==============================================================================
loginfo("Plotting figures")
gg <- (
    ggplot(data = shared_enhns)
    + geom_point(aes(x = both_DGE, y = N_shared_enhancers))
    + labs(x = "Gene pair is both differentially expressed", y = "Shared Enhancers")
    + theme_minimal()
)
ggsave("Plots/shared-enhancers.png", width = 12, height = 8, units = "cm")

# ==============================================================================
# Save data
# ==============================================================================
loginfo("Saving data")
