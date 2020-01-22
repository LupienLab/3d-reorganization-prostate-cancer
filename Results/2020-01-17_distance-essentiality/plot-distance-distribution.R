# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))

# ==============================================================================
# Constants
# ==============================================================================
metadata = fread("../../Data/External/LowC_Samples_Data_Available.tsv")
SAMPLES = paste0("PCa", metadata[, get("Sample ID")])

# ==============================================================================
# Data
# ==============================================================================
# load TSS distances for each patient
distances = rbindlist(lapply(
    SAMPLES,
    function(s) {
        dt = fread(paste0("Closest/", s, ".distance-dependency.tsv"))
        dt[, SampleID := s]
        return(dt)
    }
))

# ==============================================================================
# Analysis
# ==============================================================================
# calculate mean fraction for each gene
agg_distances = data.table(
    name = distances[, name, by = "name"]$name,
    mean_fraction = distances[, mean(Fraction), by = "name"]$V1,
    median_fraction = distances[, median(Fraction), by = "name"]$V1,
    sd_fraction = distances[, sd(Fraction), by = "name"]$V1,
    N = distances[, .N, by = "name"]$N
)

# ==============================================================================
# Plots
# ==============================================================================
gg = (
    ggplot(data = agg_distances)
    + geom_density(aes(x = mean_fraction))
    + labs(y = "Density")
    + scale_x_continuous(
        limits = c(0, 0.5),
        breaks = 0:5 / 10,
        name = "TSS distance from boundary",
        labels = paste(0:5 * 10, "%")
    )
    + theme_minimal()
)
ggsave(
    "Plots/distance-density.png",
    height = 12,
    width = 20,
    units = "cm"
)

gg = (
    ggplot(data = agg_distances)
    + stat_ecdf(aes(x = mean_fraction), geom = "step")
    + labs(y = "Cumulative Density")
    + scale_x_continuous(
        limits = c(0, 0.5),
        breaks = 0:5 / 10,
        name = "TSS distance from boundary",
        labels = paste(0:5 * 10, "%")
    )
    + theme_minimal()
)
ggsave(
    "Plots/distance-ecdf.png",
    height = 12,
    width = 20,
    units = "cm"
)
