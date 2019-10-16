# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))

# ==============================================================================
# Data
# ==============================================================================
tad_stats = fread("TADs/tad-call-stats.tsv", sep = "\t", header = TRUE)

# ==============================================================================
# Plots
# ==============================================================================
# TAD counts per patient
gg = (
    ggplot(data = tad_stats)
    + geom_point(aes(x = w, y = N, colour = Sample))
    + labs(x = "Window Size", y = bquote("Number of TADs"))
    + guides(colour = FALSE)
    + facet_grid(Sample ~ .)
)
ggsave(
    "Plots/tad-counts.png",
    height = 40,
    width = 12,
    units = "cm"
)

# TAD counts across patients
gg = (
    ggplot(data = tad_stats[, .(mean(N), sd(N)), by = w])
    + geom_errorbar(aes(x = w, ymin = V1 - V2, ymax = V1 + V2))
    + geom_point(aes(x = w, y = V1))
    + labs(x = "Window Size", y = bquote("Number of TADs"))
    + guides(colour = FALSE)
)
ggsave(
    "Plots/tad-counts-all-patients.png",
    height = 12,
    width = 20,
    units = "cm"
)

# TAD sizes per patient
gg = (
    ggplot(data = tad_stats)
    + geom_errorbar(aes(
        x = w,
        ymin = (Mean_Size - SD_Size) / 10^6,
        ymax = (Mean_Size + SD_Size) / 10^6,
        colour = Sample
    ))
    + geom_point(aes(x = w, y = Mean_Size / 10^6, colour = Sample))
    + labs(x = "Window Size", y = "Size of TADs (Mbp)")
    + guides(colour = FALSE)
    + facet_grid(Sample ~ .)
)
ggsave(
    "Plots/tad-sizes.png",
    height = 40,
    width = 12,
    units = "cm"
)
