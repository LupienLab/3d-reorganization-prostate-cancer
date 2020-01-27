# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))

# ==============================================================================
# Data
# ==============================================================================
# load read pair and dedupliation data
dedup_distance = fread("Reports/hicup_dedup_plot.tsv")
colnames(dedup_distance) = c("SampleID", "Cis_Near", "Cis_Far", "Trans", "Duplicate")

# ==============================================================================
# Analysis
# ==============================================================================
# calculate percentages based on unique read pairs
dedup_distance[, Total_Unique := Cis_Far + Cis_Far + Trans]

# convert to long form, only keep unique reads
dedup_long = melt(
    dedup_distance,
    id.vars = "SampleID",
    measure.vars = c("Cis_Near", "Cis_Far", "Trans", "Total_Unique"),
    variable.name = "Class",
    value.name = "Count"
)

# ==============================================================================
# Plots
# ==============================================================================
gg = (
    ggplot()
    + geom_col(
        data = dedup_long[Class != "Total_Unique"],
        mapping = aes(x = SampleID, y = Count, fill = Class),
        position = "stack"
    )
    + geom_text(
        data = dedup_distance,
        mapping = aes(
            x = SampleID, y = Trans / 2,
            label = paste0(round(Trans / Total_Unique, 3) * 100, "%")
        ),
        angle = 90
    )
    + labs(x = "Patient", y = "Count")
    + guides(fill = guide_legend(title = "Read pair classification"))
    + theme_minimal()
    + theme(
        axis.text.x = element_text(angle = 90),
        legend.position = "bottom"
    )
)
ggsave(
    "Reports/hicup_class.png",
    height = 12,
    width = 20,
    units = "cm"
)
