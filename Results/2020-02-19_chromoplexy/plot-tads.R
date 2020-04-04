# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))

# ==============================================================================
# Functions
# ==============================================================================

# ==============================================================================
# Data
# ==============================================================================
# load metadata
metadata <- fread(
    "../../Data/External/LowC_Samples_Data_Available.tsv",
    sep = "\t",
    header = TRUE
)

# load differential TAD results
altered_tads <- fread("Graphs/sv-disruption-tests.TADs.tsv", sep = "\t", header = TRUE)

# ==============================================================================
# Analysis
# ==============================================================================

# ==============================================================================
# Plots
# ==============================================================================
gg = (
    ggplot(data = altered_tads)
    + geom_bar(aes(x = altered_TAD, fill = altered_TAD), colour = "black")
    + labs(x = "Boundaries affected by SVs", y = "Count")
    + guides(fill = FALSE)
    + scale_x_discrete(
        breaks = c(FALSE, TRUE),
        labels = c("Not affected", "Affected")
    )
    + theme_minimal()
    + theme(
        axis.text.x = element_text(angle = 90)
    )
)
ggsave(
    "Plots/sv-disruption/tads.count.png",
    height = 12,
    width = 12,
    units = "cm",
    dpi = 400
)
ggsave(
    "Plots/sv-disruption/tads.count.pdf",
    height = 12,
    width = 12,
    units = "cm",
    dpi = 400
)
