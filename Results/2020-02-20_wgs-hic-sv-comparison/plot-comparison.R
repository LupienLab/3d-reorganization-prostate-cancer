# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
source("../2020-02-19_chromoplexy/plotting-helper.R")

# ==============================================================================
# Functions
# ==============================================================================

# ==============================================================================
# Data
# ==============================================================================
detections <- fread("detections.all.tsv", sep = "\t", header = TRUE)

# ==============================================================================
# Analysis
# ==============================================================================
detections[, Not_Mutually_Detected := Detected_In_Source - Mutually_Detected]
detections_melted <- melt(
    detections[, .SD, .SDcols = c(1, 3, 4)],
    id.vars = "Source",
    variable.name = "Detected",
    value.name = "N"
)

# ==============================================================================
# Plots
# ==============================================================================
gg = (
    ggplot(data = detections_melted)
    + geom_col(
        aes(x = Source, y = N, fill = Detected),
        position = position_stack(),
        colour = "#000000"
    )
    + geom_text(
        aes(x = Source, y = N, label = paste0(N, " (", 100 * round(N / detections[, Detected_In_Source], 3), "%)"), group = Detected),
        position = position_stack(vjust = 0.5)
    )
    + labs(x = "Detection In", y = "Breakpoints Detected")
    + scale_fill_manual(
        name = "Detected in Other",
        limits = c("Not_Mutually_Detected", "Mutually_Detected"),
        labels = c("No", "Yes"),
        values = c("#bdbdbd", "#f8766d")
    )
    + guides(colour = FALSE)
    + theme_minimal()
    + theme(
        legend.position = "bottom"
    )
)
savefig(gg, "detections.all")
