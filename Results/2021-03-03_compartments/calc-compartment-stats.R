# ==============================================================================
# Meta
# ==============================================================================
# calc-compartment-stats
# ------------------------------------------------
# Author: James Hawley
# Description: Calculate the statistics of compartment calls for a single sample


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

# sample metadata
meta <- fread("config.tsv", sep = "\t", header = TRUE)
meta <- meta[Include == "Yes"]
SAMPLES <- meta[, Sample_ID]

compartments <- rbindlist(lapply(
    SAMPLES,
    function(s) {
        dt <- fread(
            file.path(CMPMT_DIR, paste0(s, ".compartments.cis.vecs.tsv")),
            sep = "\t",
            header = TRUE
        )
        dt[, Sample_ID := s]
        return(dt)
    }
))

# ==============================================================================
# Analysis
# ==============================================================================
loginfo("Calculating statistics")

# 1. identify the compartment boundaries by detecting a sign switch in E1 between consecutive rows
# ------------------------------------------------
# first remove regions with NA values
compartments <- compartments[complete.cases(compartments)]

# find bandwidth for appropriate determination of compartment boundaries
bandwidth <- compartments[,
    .(diff = E1 - shift(E1, 1)),
    by = c("Sample_ID", "chrom")
]
cmpmt_diff_threshold <- 1.5 * bandwidth[, sd(diff, na.rm = TRUE)]

compartments <- compartments[,
    .(
        start,
        end,
        E1,
        new_compartment_start = (
            (sign(E1) != sign(shift(E1, 1)))
            & (abs(E1 - shift(E1, 1)) > cmpmt_diff_threshold)
        )
    ),
    by = c("Sample_ID", "chrom")
]
# NAs means the start of a new chromosome and should be changed to TRUE
compartments[is.na(new_compartment_start), new_compartment_start := TRUE]

# calculate number of compartments per sample
cmpmt_counts <- compartments[new_compartment_start == TRUE, .N, by = "Sample_ID"]
cmpmt_counts[, summary(N)]

# calculate the size of compartments
cmpmt_sizes <- compartments[
    new_compartment_start == TRUE,
    .(size = start - shift(start, 1)),
    by = c("Sample_ID", "chrom")
]
# remove the NAs, which are the start of new chromosomes or samples
cmpmt_sizes <- cmpmt_sizes[complete.cases(cmpmt_sizes)]

# merge with sample metadata
cmpmt_sizes <- merge(
    x = cmpmt_sizes,
    y = meta[, .SD, .SDcols = c("Sample_ID", "Type")],
    by = "Sample_ID"
)
cmpmt_sizes[, mean(size)]

# ==============================================================================
# Plots
# ==============================================================================
loginfo("Plotting figures")

gg_n <- (
    ggplot(
        data = cmpmt_counts,
        mapping = aes(
            x = Sample_ID,
            y = N,
            fill = Sample_ID
        )
    )
    + geom_col()
    + scale_x_discrete(
        name = NULL,
        breaks = meta[, Sample_ID],
        labels = meta[, Label]
    )
    + scale_y_continuous(
        name = "Compartments"
    )
    + scale_fill_manual(
        breaks = meta[, Sample_ID],
        labels = meta[, Label],
        values = meta[, Type_Colour]
    )
    + guides(fill = FALSE)
    + coord_flip()
    + theme_minimal()
    + theme(
        axis.text.x = element_text(colour = "#000000"),
        axis.text.y = element_text(colour = "#000000"),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "#000000"),
        panel.grid.minor.x = element_line(colour = "#000000")
    )
)
ggsave(
    file.path(PLOT_DIR, "compartment-counts.png"),
    gg_n,
    width = 11,
    height = 12,
    units = "cm"
)
ggsave(
    file.path(PLOT_DIR, "compartment-counts.pdf"),
    gg_n,
    width = 11,
    height = 12,
    units = "cm"
)

gg_size <- (
    ggplot(
        data = cmpmt_sizes,
        mapping = aes(
            x = size / 1000,
            fill = Sample_ID
        )
    )
    + geom_density(aes(y = ..count..), alpha = 0.1)
    + scale_x_continuous(
        name = "Compartment Size (kbp)"
        # limits = c(0, 1e3)
    )
    + scale_y_continuous(
        name = "Density"
    )
    + scale_fill_manual(
        breaks = meta[, Sample_ID],
        labels = meta[, Label],
        values = meta[, Type_Colour]
    )
    + facet_grid(. ~ Type)
    + guides(fill = FALSE)
    + theme_minimal()
    + theme(
        axis.text.x = element_text(colour = "#000000"),
        axis.text.y = element_text(colour = "#000000")
    )
)
ggsave(
    file.path(PLOT_DIR, "compartment-sizes.png"),
    gg_size,
    width = 12,
    height = 8,
    units = "cm"
)
ggsave(
    file.path(PLOT_DIR, "compartment-sizes.pdf"),
    gg_size,
    width = 12,
    height = 8,
    units = "cm"
)

# ==============================================================================
# Save Data
# ==============================================================================
loginfo("Saving data")
fwrite(
    cmpmt_counts,
    file.path(CMPMT_DIR, "compartment-counts.tsv"),
    sep = "\t",
    col.names = TRUE
)
fwrite(
    cmpmt_sizes[, .N, by = c("Sample_ID", "size")],
    file.path(CMPMT_DIR, "compartment-sizes.tsv"),
    sep = "\t",
    col.names = TRUE
)

