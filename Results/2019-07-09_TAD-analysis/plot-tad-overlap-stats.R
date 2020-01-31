# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))

# ==============================================================================
# Data
# ==============================================================================
# read overlapping statistics
samplewise_intersections = fread("TAD-comparisons/comparison-total-counts.tsv", sep = "\t", header = TRUE)
WINDOWS = 3:30

# ==============================================================================
# Plots
# ==============================================================================
# convert to factor for ordered plotting
samplewise_intersections[, Window := factor(Window, levels = WINDOWS, ordered = TRUE)]

gg = (
    ggplot(data = samplewise_intersections)
    + geom_boxplot(aes(x = Window, y = 100 * Frac, fill = Window, group = Window), outlier.shape = NA)
    + labs(x = "Window Size", y = "Similar TADs (%)")
    + guides(colour = FALSE, fill = FALSE)
    + scale_x_discrete(
        breaks = c(3, 10, 20, 30),
        labels = c(3, 10, 20, 30),
        limits = c(3, 10, 20, 30)
    )
    + scale_fill_viridis_c()
    + facet_wrap(~ Sample1, ncol = 7)
    + theme_minimal()
    + theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        # panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        # panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()
    )
)
ggsave(
    "Plots/tad-similarity-counts.png",
    height = 12,
    width = 30,
    units = "cm"
)
