# ==============================================================================
# Meta
# ==============================================================================
# plot-compartments
# ------------------------------------------------
# Author: James Hawley
# Description: Plot various plots about compartments


# ==============================================================================
# Environment
# ==============================================================================
suppressWarnings(library("logging"))
loginfo("Loading packages")

suppressWarnings(library("data.table"))
suppressWarnings(library("ggplot2"))


# ==============================================================================
# Functions
# ==============================================================================


# ==============================================================================
# Data
# ==============================================================================
loginfo("Loading data")


# ==============================================================================
# Plots
# ==============================================================================
loginfo("Plotting figures")

gg <- (
    ggplot(data = wide_e1)
    + geom_point(
        mapping = aes(x = All_Mean, y = All_SD ^ 2),
        alpha = 0.01
    )
    + scale_x_continuous(
        name = "Mean"
    )
    + scale_y_continuous(
        name = "Variance"
    )
    + theme_minimal()
)
ggsave(
    "e1-mean-sd.png",
    gg,
    width = 8,
    height = 8,
    units = "cm"
)

gg <- (
    ggplot(data = wide_e1)
    + geom_point(
        mapping = aes(
            x = (Benign_Mean + Malignant_Mean) / 2,
            y = (Benign_Mean - Malignant_Mean)
        ),
        alpha = 0.01
    )
    + scale_x_continuous(
        name = "Mean"
    )
    + scale_y_continuous(
        name = "Difference"
    )
    + theme_minimal()
)
ggsave(
    "e1-mean-diff.png",
    gg,
    width = 8,
    height = 8,
    units = "cm"
)

gg <- (
    ggplot(data = wide_e1)
    + geom_abline(
        intercept = 0,
        slope = 1,
        linetype = "dashed"
    )
    + geom_point(
        mapping = aes(
            x = (Benign_Mean + Malignant_Mean) / 2,
            y = Malignant_Mean, size = Benign_SD),
        alpha = 0.1
    )
    + theme_minimal()
)
ggsave(
    "e1.png",
    gg,
    width = 12,
    height = 12,
    units = "cm"
)



# ==============================================================================
# Save Data
# ==============================================================================
loginfo("Saving data")

