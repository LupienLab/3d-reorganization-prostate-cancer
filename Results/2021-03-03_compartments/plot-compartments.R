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


CMPMT_DIR = "Compartments"
PLOT_DIR = "Plots"

# ==============================================================================
# Data
# ==============================================================================
loginfo("Loading data")
wide_e1 <- fread(
    file.path(CMPMT_DIR, "compartments.stats.tsv"),
    sep = "\t",
    header = TRUE
)

# ==============================================================================
# Plots
# ==============================================================================
loginfo("Plotting figures")

gg <- (
    ggplot(data = wide_e1)
    + geom_point(
        mapping = aes(x = All_Mean, y = All_SD ^ 2),
        alpha = 0.01,
        shape = 21
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
    file.path(PLOT_DIR, "e1-mean-sd.png"),
    gg,
    width = 8,
    height = 8,
    units = "cm"
)
ggsave(
    file.path(PLOT_DIR, "e1-mean-sd.pdf"),
    gg,
    width = 8,
    height = 8,
    units = "cm"
)

gg_chr3 <- (
    ggplot(data = wide_e1[chrom == "chr3"])
    + geom_point(
        mapping = aes(
            x = (Benign_Mean + Malignant_Mean) / 2,
            y = (Benign_Mean - Malignant_Mean)
        ),
        alpha = 0.01,
        shape = 21
    )
    + scale_x_continuous(
        name = "(Benign + Tumour) / 2",
        limits = c(-2, 2)
    )
    + scale_y_continuous(
        name = "Benign - Tumour",
        limits = c(-2, 1)
    )
    + theme_minimal()
)
ggsave(
    file.path(PLOT_DIR, "e1-mean-diff.chr3.png"),
    gg,
    width = 8,
    height = 8,
    units = "cm"
)
ggsave(
    file.path(PLOT_DIR, "e1-mean-diff.chr3.pdf"),
    gg,
    width = 8,
    height = 8,
    units = "cm"
)

gg_chr19 <- (
    ggplot(data = wide_e1[chrom == "chr19"])
    + geom_point(
        mapping = aes(
            x = (Benign_Mean + Malignant_Mean) / 2,
            y = (Benign_Mean - Malignant_Mean)
        ),
        alpha = 0.01,
        shape = 21
    )
    + scale_x_continuous(
        name = "(Benign + Tumour) / 2",
        limits = c(-2, 2)
    )
    + scale_y_continuous(
        name = "Benign - Tumour",
        limits = c(-2, 1)
    )
    + theme_minimal()
)
ggsave(
    file.path(PLOT_DIR, "e1-mean-diff.chr19.png"),
    gg,
    width = 8,
    height = 8,
    units = "cm"
)
ggsave(
    file.path(PLOT_DIR, "e1-mean-diff.chr19.pdf"),
    gg,
    width = 8,
    height = 8,
    units = "cm"
)

gg_all <- (
    ggplot(
        data = wide_e1,
        mapping = aes(
            x = (Benign_Mean + Malignant_Mean) / 2,
            y = (Benign_Mean - Malignant_Mean),
            fill = chrom,
            colour = chrom
        )
    )
    + geom_point(
        data = subset(wide_e1, ),
        alpha = 0.1,
        shape = 21,
        size = 2
    )
    + geom_point(
        data = subset(wide_e1, (chrom == "chr3")),
        alpha = 0.1,
        shape = 21,
        size = 2
    )
    + geom_point(
        # data = subset(wide_e1, (chrom == "chr19") & (start > 20000000)),
        data = subset(wide_e1, chrom == "chr19"),
        alpha = 0.2,
        shape = 21,
        size = 2
    )
    + geom_point(
        # data = subset(wide_e1, (chrom == "chrY") & (start > 10000000) & (end < 16000000)),
        data = subset(wide_e1, chrom == "chrY"),
        alpha = 0.2,
        shape = 21,
        size = 2
    )
    + scale_x_continuous(
        name = "Compartmentalization Mean",
        limits = c(-2, 2)
    )
    + scale_y_continuous(
        name = "Compartmentalization Difference (Benign - Tumour)",
        limits = c(-2, 1)
    )
    + scale_colour_manual(
        breaks = paste0("chr", c(1:22, "X", "Y")),
        values = c(
            rep("#000000", 2),
            "#00A9FF",
            rep("#000000", 15),
            "#00C08B",
            rep("#000000", 4),
            "#FF6C91"
        )
    )
    + scale_fill_manual(
        breaks = paste0("chr", c(1:22, "X", "Y")),
        values = c(
            rep("#000000", 2),
            "#00A9FF",
            rep("#000000", 15),
            "#00C08B",
            rep("#000000", 4),
            "#FF6C91"
        )
    )
    + guides(
        colour = FALSE,
        fill = guide_legend(override.aes = list(alpha = 1))
    )
    + theme_minimal()
    + theme(
        legend.position = "bottom",
        axis.text = element_text(colour = "#000000"),
    )
)
ggsave(
    file.path(PLOT_DIR, "e1-mean-diff.png"),
    gg_all,
    width = 20,
    height = 20,
    units = "cm"
)
ggsave(
    file.path(PLOT_DIR, "e1-mean-diff.pdf"),
    gg_all,
    width = 20,
    height = 20,
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
            y = Malignant_Mean, size = Benign_SD
        ),
        alpha = 0.1,
        shape = 21
    )
    + theme_minimal()
)
ggsave(
    file.path(PLOT_DIR, "e1.png"),
    gg,
    width = 12,
    height = 12,
    units = "cm"
)
ggsave(
    file.path(PLOT_DIR, "e1.pdf"),
    gg,
    width = 12,
    height = 12,
    units = "cm"
)
