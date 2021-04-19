# ==============================================================================
# Meta
# ==============================================================================
# Plot loops
# --------------------------------------
# Description: Plot statistical information about loop calls
# Author: James Hawley


# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("logging"))

loginfo("Loading packages")
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))

LOOP_DIR <- "Loops"
PLOT_DIR <- "Plots"

# ==============================================================================
# Data
# ==============================================================================
loginfo("Loading data")

# load metadata
meta <- fread("config.tsv", sep = "\t", header = TRUE)
meta <- meta[Include == "Yes"]

# load loop calls
loops_indiv <- fread(
    file.path(LOOP_DIR, "merged-loops.tsv"),
    sep = "\t",
    header = TRUE
)
loops <- fread(
    file.path(LOOP_DIR, "merged-loops.sample-counts.tsv"),
    sep = "\t",
    header = TRUE
)

# ==============================================================================
# Analysis
# ==============================================================================
loginfo("Performing calculations")

# merge loop counts and sample metadata
loop_counts <- merge(
    x = loops_indiv[, .N, by = "SampleID"],
    y = meta[, .SD, .SDcols = c("SampleID", "Type")],
    by = "SampleID"
)

# calculate distances between loop anchors
loops[, Distance := pmin(abs(start_x - end_y), abs(end_x - start_y))]

distance_stats <- data.table(
    Statistic = c("Mean", "Median", "Mode"),
    Value = c(
        loops[, mean(Distance)],
        loops[, median(Distance)],
        loops[, as.numeric(names(which.max(table(Distance))))]
    )
)

# recast inidividual loop calls for upset plot
loops_upset <- dcast(
    loops_indiv,
    loop_ID ~ SampleID,
    value.var = "SampleID",
    fun.aggregate = function(x) length(x) > 0
)
loops_grid <- melt(
    loops_upset,
    id.vars = "loop_ID",
    variable.name = "SampleID",
    value.name = "contains"
)
# merge metadata
loops_grid <- merge(
    x = loops_grid,
    y = meta[, .SD, .SDcols = c("SampleID", "Label", "Type")]
)

# calculate the overlap of loops by type, and the number of each tissue type with this loop
loops_grid_by_type <- dcast(
    loops_grid,
    loop_ID ~ Type,
    value.var = "contains",
    fun.aggregate = function(x) length(which(x))
)
loops_grid_by_type[, `:=`(
    n_intersections = Benign + Malignant,
    type_intersection = ifelse(
        Benign > 0,
        ifelse(
            Malignant > 0,
            "Shared",
            "Benign-Specific"
        ),
        "Tumour-Specific"
    )
)]
loops_grid_by_type_summary <- loops_grid_by_type[, .N, by = c("type_intersection", "n_intersections")]

# ==============================================================================
# Plots
# ==============================================================================
loginfo("Plotting data")
gg <- (
    ggplot(data = loops)
    + geom_density(aes(x = Distance / 1000))
    + scale_x_continuous(
        name = "Distance between loop anchors (kbp)"
    )
    + geom_vline(
        data = distance_stats,
        mapping = aes(xintercept = Value / 1000, colour = Statistic),
        linetype = "dashed"
    )
    + geom_text(
        data = distance_stats,
        mapping = aes(
            x = Value / 1000,
            y = 0,
            colour = Statistic,
            label = paste(Statistic, "=", round(Value / 1000))
        ),
        angle = 90,
        hjust = 0,
        vjust = 1.1
    )
    + scale_colour_manual(
        limits = c("Mean", "Median", "Mode"),
        values = c("#1e90ff", "#66cdaa", "#b22222")
    )
    + guides(colour = FALSE)
    + theme_minimal()
)
ggsave("Plots/loop-calls.anchor-distance.png", gg, width = 12, height = 8, units = "cm")

gg_loops_samples <- (
    ggplot(
        data = loop_counts,
        mapping = aes(
            x = SampleID,
            y = N,
            fill = SampleID
        )
    )
    + geom_col()
    + scale_x_discrete(
        name = NULL,
        breaks = meta[, SampleID],
        labels = meta[, Label]
    )
    + scale_y_continuous(
        name = "Interactions Called"
    )
    + scale_fill_manual(
        breaks = meta[, SampleID],
        values = meta[, Type_Colour]
    )
    + guides(fill = FALSE)
    + coord_flip()
    + theme_minimal()
    + theme(
        axis.text = element_text(colour = "#000000"),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "#000000"),
        panel.grid.minor.x = element_line(colour = "#000000")
    )
)
ggsave(
    file.path("Plots", "loop-calls.per-sample.png"),
    gg_loops_samples,
    width = 12,
    height = 8,
    units = "cm"
)
ggsave(
    file.path("Plots", "loop-calls.per-sample.pdf"),
    gg_loops_samples,
    width = 12,
    height = 8,
    units = "cm"
)


gg_calls_by_type_all <- (
    ggplot(
        data = loops_grid_by_type_summary,
        mapping = aes(x = n_intersections, y = N, fill = type_intersection)
    )
    + scale_x_continuous(
        name = "Number of Samples with Interaction"
    )
    + scale_y_continuous(
        name = "Frequency"
    )
    + scale_fill_manual(
        name = "Interaction Type",
        breaks = c(
            "Benign-Specific",
            "Shared",
            "Tumour-Specific"
        ),
        values = c(
            "#AEC7E8",
            "#B5B5B5",
            "#2077B4"
        )
    )
    + geom_col(position = position_stack())
    + theme_minimal()
    + theme(
        legend.position = "bottom",
        axis.text = element_text(colour = "#000000")
    )
)
ggsave(
    file.path(PLOT_DIR, "loop-calls.by-type.all.png"),
    gg_calls_by_type_all,
    width = 20,
    height = 8,
    units = "cm"
)
ggsave(
    file.path(PLOT_DIR, "loop-calls.by-type.all.pdf"),
    gg_calls_by_type_all,
    width = 20,
    height = 8,
    units = "cm"
)
