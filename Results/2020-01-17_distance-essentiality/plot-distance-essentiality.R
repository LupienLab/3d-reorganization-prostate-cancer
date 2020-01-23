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
# load distance calculations for each patient
distances = rbindlist(lapply(
    SAMPLES,
    function(s) {
        dt = fread(paste0("Closest/", s, ".distance-dependency.tsv"))
        dt[, SampleID := s]
        return(dt)
    }
))

# load gene essentiality information for a few prostate cell lines
depmap = fread("../../Data/External/DepMap/depmap-rnai.tsv")
# melt into long form
depmap_long = melt(
    depmap,
    id.vars = "Gene",
    variable.name = "Cell",
    value.name = "Essentiality"
)

cells = colnames(depmap)[-1]

# ==============================================================================
# Analysis
# ==============================================================================
# calculate median and mean distances over patients
agg_distances = distances[, .(mean(Fraction), sd(Fraction)), by = "name"]
colnames(agg_distances) = c("Gene", "Mean_Fraction", "SD_Fraction")

# merge depmap and distance information
# the DepMap data contains 17 309 genes, and the distance data from the GENCODE annotation contains 19 926
# merging the two by their HUGO gene name leads to 15 922
combined = merge(
    x = agg_distances,
    y = depmap_long,
    by = "Gene"
)

# calculate the quantile each gene belongs to, wrt the individual
for (cell in cells) {
    breaks = depmap[, quantile(get(cell), 0:5 / 5, na.rm = TRUE)]
    combined[
        Cell == cell,
        Level := cut(
            Essentiality,
            breaks = breaks,
            labels = as.character(1:5)
        )
    ]
}

# ==============================================================================
# Plots
# ==============================================================================
# plot essentiality across TADs
gg = (
    ggplot(data = combined[complete.cases(combined)])
    + geom_density(aes(x = Mean_Fraction, colour = Level))
    + labs(y = "Density")
    + scale_x_continuous(
        limits = c(0, 0.5),
        breaks = 0:5 / 10,
        name = "Mean TSS distance from boundary",
        labels = paste0(0:5 * 10, "%")
    )
    + scale_colour_manual(
        breaks = 1:5,
        name = "Essentiality Quantile",
        # reversing these lists because according to DepMap
        # the more negative a number, the more essential that gene
        labels = rev(c(
            "[0%, 20%)",
            "[20%, 40%)",
            "[40%, 60%)",
            "[60%, 80%)",
            "[80%, 100%]"
        )),
        values = rev(c(
            "#fef0d9",
            "#fdcc8a",
            "#fc8d59",
            "#e34a33",
            "#b30000"
        ))
    )
    + facet_wrap(~ Cell)
    + theme_minimal()
    + theme(
        axis.text.x = element_text(angle = 90)
    )
)
ggsave(
    "Plots/distance-density-by-essentiality.png",
    height = 12,
    width = 20,
    units = "cm"
)
