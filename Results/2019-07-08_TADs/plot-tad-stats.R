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
gg = (
    ggplot(data = tad_stats)
    + geom_point(aes(x = w, y = N, colour = Sample))
    + labs(x = "Window Size", y = bquote("Number of TADs"))
    + guides(colour = FALSE)
    + facet_wrap(~ Sample)
)
ggsave(
    "Plots/tad-counts.png",
    height = 12,
    width = 20,
    units = "cm"
)

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
    + facet_wrap(~ Sample)
)
ggsave(
    "Plots/tad-sizes.png",
    height = 12,
    width = 20,
    units = "cm"
)
