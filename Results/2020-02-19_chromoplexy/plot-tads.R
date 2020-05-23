# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
source("plotting-helper.R")

# ==============================================================================
# Data
# ==============================================================================
# load metadata
metadata <- fread(
    "../../Data/External/LowC_Samples_Data_Available.tsv",
    sep = "\t",
    header = TRUE
)
metadata <- metadata[Include == "Yes"]

# load differential TAD results
altered_tads <- fread("Graphs/sv-disruption-tests.TADs.tsv", sep = "\t", header = TRUE)

# ==============================================================================
# Plots
# ==============================================================================
gg = (
    ggplot()
    + geom_bar(
        data = altered_tads,
        aes(x = altered_TAD, fill = altered_TAD),
        colour = "black"
    )
    + geom_text(
        data = altered_tads[, .N, by = "altered_TAD"],
        aes(x = altered_TAD, y = N, label = paste0(N, " (", 100 * round(N / sum(N), 3), "%)")),
        colour = "black",
        vjust = -0.5
    )
    + labs(x = "SV affects TAD boundaries", y = "Count")
    + guides(fill = FALSE)
    + scale_x_discrete(
        breaks = c(FALSE, TRUE),
        labels = c("No", "Yes")
    )
    + scale_fill_manual(
        breaks = c(FALSE, TRUE),
        labels = c("No", "Yes"),
        values = c("#BDBDBD", "#000000")
    )
    + theme_minimal()
    + theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    )
)
savefig(gg, "Plots/sv-disruption/tads.count", width = 12)
