# Plot the distance for each TAD boundary to a CTCF motif
#
# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))

# ==============================================================================
# Data
# ==============================================================================
# read intersection data
dt <- fread(
    "Proximal/nearest-motif.bed",
    sep = "\t",
    header = FALSE,
    col.names = c(
        "chr_bound", "start_bound", "end_bound", "N_Int", "Which_Int",
        1:13,
        "chr_motif", "start_motif", "end_motif", "strand_motif", "motifs",
        "distance"
    )
)

# ==============================================================================
# Analysis
# ==============================================================================

# ==============================================================================
# Plots
# ==============================================================================
gg <- (
    ggplot(data = dt)
    # + 1 to not discard boundaries that overlap with motifs
    + geom_density(aes(x = distance + 1))
    + labs(
        x = "Distance between CTCF motif and TAD boundary (bp)",
        y = "Density"
    )
    + scale_x_log10()
    + theme_minimal()
)
ggsave(
    "Plots/boundary-motif-distance.png",
    height = 12,
    width = 20,
    units = "cm"
)

# stratify densities by how frequently they are shared between patients
gg <- (
    ggplot(data = dt)
    # + 1 to not discard boundaries that overlap with motifs
    + geom_density(aes(x = distance + 1, fill = factor(N_Int)), alpha = 0.3)
    + labs(
        x = "Distance between CTCF motif and TAD boundary (bp)",
        y = "Density"
    )
    + scale_x_log10()
    + guides(fill = guide_legend(title = "Shared by"))
    + theme_minimal()
)
ggsave(
    "Plots/boundary-motif-distance.stratified.png",
    height = 12,
    width = 20,
    units = "cm"
)
