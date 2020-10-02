# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
source("../2020-02-19_chromoplexy/plotting-helper.R")

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
altered_tads <- fread("sv-disruption-tests.TADs.tsv", sep = "\t", header = TRUE)

altered_tads_count <- altered_tads[, .N, by = "altered_TAD"]

# ==============================================================================
# Plots
# ==============================================================================
gg <- (
    ggplot(data = altered_tads_count)
    + geom_col(
        aes(x = 1, y = N, fill = factor(altered_TAD, ordered = TRUE, levels = c(TRUE, FALSE))),
        colour = "black",
        position = "stack"
    )
    + geom_text(
        aes(
            x = 1,
            y = N,
            label = paste0(N, " (", 100 * round(N / sum(N), 3), "%)")
        ),
        colour = "white",
        position = position_stack(vjust = 0.5)
    )
    + labs(x = NULL, y = "Breakpoints affecting\nTAD boundaries")
    + guides(fill = FALSE)
    + scale_fill_manual(
        breaks = c(TRUE, FALSE),
        labels = c("Yes", "No"),
        values = c("#000000", "#AEA28E")
    )
    + theme_minimal()
    + theme(
        legend.position = "bottom",
        axis.text.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank()
    )
)
savefig(gg, "Plots/tads.count", width = 5, height = 5.5)

# ==============================================================================
# Save data
# ==============================================================================
fwrite(
    altered_tads_count,
    "sv-disruption-tests.TADs.counts.tsv",
    sep = "\t",
    col.names = TRUE
)
 