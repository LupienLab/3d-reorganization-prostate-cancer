# ==============================================================================
# Meta
# ==============================================================================
# sv-complexity
# --------------------------------------
# Description: Assess whether the SV's functionality is related to its size and complexity
# Author: James Hawley


# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("logging"))

loginfo("Loading packages")
suppressMessages(library("data.table"))
suppressMessages(library("GenomicRanges"))
suppressMessages(library("ggplot2"))
source(file.path("..", "src", "savefig.R"))


# ==============================================================================
# Data
# ==============================================================================
loginfo("Loading data")

# load metadata
metadata <- fread("config.tsv", sep = "\t")

# load SV breakpoints and their expression changes
sv_n_loops <- fread("sv-loop-intersection.tsv", sep = "\t")
sv_components <- fread("../2020-02-19_chromoplexy/Statistics/sv-components.tsv", sep = "\t")
sv_exprs <- fread("../2020-06-18_sv-disruption-expression/1sample-results.tsv", sep = "\t")

# load loop anchors
loops <- fread("../2020-10-06_loops/Loops/merged-loops.sample-counts.tsv", sep = "\t")
anchors <- rbind(
    loops[, .(
        chr = chr_x,
        start = start_x,
        end = end_x,
        anchor_ID = anchor_ID_x,
        loop_ID = loop_ID,
        fdr = fdr,
        detection_scale = detection_scale
    )],
    loops[, .(
        chr = chr_y,
        start = start_y,
        end = end_y,
        anchor_ID = anchor_ID_y,
        loop_ID = loop_ID,
        fdr = fdr,
        detection_scale = detection_scale
    )]
)

# ==============================================================================
# Analysis
# ==============================================================================
loginfo("Calculating observed intersection of SV breakpoints and loops")

# merge SV complexity into previous table
sv_n_loops <- merge(
    x = sv_n_loops,
    y = sv_components[, .SD, .SDcols = c("SampleID", "component_ID", "N_Breakpoints", "N_Chr")],
    by = c("SampleID", "component_ID")
)

# test for correlation between number of loops and complexity of SV
testing_cor <- cor.test(
    x = sv_n_loops$n_loops,
    y = sv_n_loops$N_Breakpoints,
    method = "spearman",
    alternative = "greater"
)

# ==============================================================================
# Plots
# ==============================================================================
loginfo("Making figures")
gg <- (
    ggplot(data = sv_n_loops)
    + geom_point(aes(x = N_Breakpoints, y = n_loops))
    + labs(x = "Breakpoints", y = "Loops")
    + theme_minimal()
)
savefig(gg, "Plots/loops-sv-complexity")
