# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
source(
    file.path("..", "2020-02-19_chromoplexy", "plotting-helper.R")
)

RES_DIR <- file.path("..", "..", "results", "2020-02-20_wgs-hic-sv-comparison")
PLOT_DIR <- file.path(RES_DIR, "Plots")

# ==============================================================================
# Data
# ==============================================================================
# laod metadata
meta <- fread(
    file.path("..", "config.tsv"),
    sep = "\t"
)
# retain primary tumour samples
meta <- meta[(Source == "Primary") & (Type == "Malignant")]

detections <- fread(
    file.path(RES_DIR, "detections.all.tsv"),
    sep = "\t",
    header = TRUE
)
detections_sample <- fread(
    file.path(RES_DIR, "detections.per-sample.tsv"),
    sep = "\t",
    header = TRUE
)
sv_types <- fread(    
    file.path(RES_DIR, "sv-types-counted.tsv"),
    sep = "\t",
    header = TRUE
)

# ==============================================================================
# Analysis
# ==============================================================================
detections[, Not_Mutually_Detected := Detected_In_Source - Mutually_Detected]
detections_sample[, Not_Mutually_Detected := Detected_In_Source - Mutually_Detected]
detections_melted <- melt(
    detections[, .SD, .SDcols = c(1, 3, 4)],
    id.vars = "Source",
    variable.name = "Detected",
    value.name = "N"
)
detections_sample_melted <- melt(
    detections_sample[, .SD, .SDcols = c(1, 2, 4, 5)],
    id.vars = c("SampleID", "Source"),
    variable.name = "Detected",
    value.name = "N"
)

# ==============================================================================
# Plots
# ==============================================================================
gg_all <- (
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
savefig(gg_all, file.path(PLOT_DIR, "detections.all"))

gg_sample <- (
    ggplot(data = detections_sample_melted)
    + geom_col(
        aes(x = Source, y = N, fill = Detected),
        position = position_stack(),
        colour = "#000000"
    )
    + labs(x = "Detection In", y = "Breakpoints Detected")
    + scale_fill_manual(
        name = "Detected in Other",
        limits = c("Not_Mutually_Detected", "Mutually_Detected"),
        labels = c("No", "Yes"),
        values = c("#bdbdbd", "#f8766d")
    )
    + guides(colour = FALSE)
    + facet_wrap(~ SampleID, nrow = 1)
    + theme_minimal()
    + theme(
        legend.position = "bottom",
        axis.text.x = element_text(angle = 90)
    )
)
savefig(
    gg_sample,
    file.path(PLOT_DIR, "detections.per-patient"),
    width = 30
)

# helper function for nicely plotting facets
map_to_label <- function(s) {
    plot_label <- meta[Sample_ID == s, Label]
    ifelse(
        # if no label, return a blank value
        nchar(plot_label) == 0,
        "",
        # otherwise return the "Label"
        plot_label
    )
}

gg_sv_types <- (
    ggplot(
        data = sv_types,
        mapping = aes(
            x = SV_Type,
            y = N,
            fill = Method
        )
    )
    + geom_col(
        position = position_dodge(),
        colour = "#000000"
    )
    + facet_wrap(
        ~ SampleID,
        nrow = 2,
        labeller = labeller(SampleID = map_to_label)
    )
    + scale_y_continuous(
        name = "Frequency"
    )
    + theme_minimal()
    + theme(
        legend.position = "bottom",
        axis.text.x = element_text(
            angle = 90,
            vjust = 0.5,
            hjust = 1,
            colour = "#000000"
        ),
        axis.ticks.x = element_blank(),
        axis.ticks.length.y = unit(4, "pt"),
        axis.ticks.y = element_line(
            colour = "#000000"
        ),
        axis.line.y.left = element_line(
            colour = "#000000"
        ),
        axis.line.x.bottom = element_blank(),
        panel.grid = element_blank()
    )
)
savefig(
    gg_sv_types,
    file.path(PLOT_DIR, "detections.typed")
)
