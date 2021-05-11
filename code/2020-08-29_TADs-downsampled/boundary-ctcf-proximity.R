# ==============================================================================
# Meta
# ==============================================================================
# Boundary CTCF proximity
# --------------------------------------
# Description: Calculate enrichment of CTCF binding sites near TAD boundaries
# Author: James Hawley


# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("logging"))

loginfo("Loading packages")
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
suppressMessages(library("argparse"))

if (!interactive()) {
    PARSER <- argparse::ArgumentParser(
        description = "Calculate enrichment of CTCF binding sites near TAD boundaries"
    )
    PARSER$add_argument(
        "cell",
        type = "character",
        help = "Cell line with CTCF binding sites to compare to",
        choices = c("LNCaP", "22Rv1", "C4-2B", "VCaP")
    )
    ARGS <- PARSER$parse_args()
} else {
    ARGS <- list(
        "cell" = "22Rv1"
    )
}


# ==============================================================================
# Data
# ==============================================================================
loginfo("Loading data")
metadata <- fread("config.tsv", sep = "\t", header = TRUE)
metadata <- metadata[Include == "Yes"]
ALL_SAMPLES <- metadata[, SampleID]

# load boundary-CTCF peak pairings for each sample
pairs <- rbindlist(lapply(
    ALL_SAMPLES,
    function(s) {
        dt <- fread(
            file.path("CTCF", "Intersection", paste0(s, ".", ARGS$cell, ".CTCF-peaks.bed")),
            sep = "\t",
            header = FALSE,
            # ignore intersection,
            select = c(1:3, 15),
            col.names = c(
                "chr_peak", "start_peak", "end_peak", "pos_bound"
            )
        )
        dt[, SampleID := s]
        return(dt)
    }
))


# ==============================================================================
# Analysis
# ==============================================================================
loginfo("Calculating frequencies")
# calculate distance between CTCF peak and the TAD boundary
pairs[, Distance := ifelse(
    (start_peak <= pos_bound) & (pos_bound <= end_peak),
    0,
    sign(start_peak - pos_bound) * pmin(
        abs(start_peak - pos_bound),
        abs(end_peak - pos_bound)
    )
)]

# create 5 kbp bins around the boundary
pairs[, Distance_Bin := cut(
    Distance,
    seq(-200000 + 2500, 200000 - 2500, by = 5000)
)]

# calculate
all_counts <- merge(
    # calculate number of peaks at a distance bin, per patient
    x = pairs[, .(N_Bin = .N), keyby = c("SampleID", "Distance_Bin")],
    # calculate number of peaks in total, per patient
    y = pairs[, .(N_Total = .N), keyby = c("SampleID")],
    by = "SampleID"
)
# calculate the fraction at each distance bin
all_counts[, Freq := N_Bin / N_Total]

# extract upper and lower bound from cut factor labels
all_counts[, Bin_Bounds := gsub("\\[|\\]|\\(|\\)", "", Distance_Bin)]
all_counts[, Bin_Lower := as.numeric(gsub(",[-e0-9\\.\\+]+", "", Bin_Bounds))]
all_counts[, Bin_Upper := as.numeric(gsub("[-e0-9\\.\\+]+,", "", Bin_Bounds))]
all_counts[, Bin_Mid := (Bin_Lower + Bin_Upper) / 2]


# ==============================================================================
# Save data
# ==============================================================================
loginfo("Saving data")
fwrite(
    all_counts[,
        .SD,
        .SDcols = c(
            "SampleID", "Distance_Bin", "N_Bin", "N_Total",
            "Freq", "Bin_Lower", "Bin_Mid", "Bin_Upper"
        )
    ],
    file.path("CTCF", "Distances", paste0(ARGS$cell, "-CTCF-peaks.distances.tsv")),
    sep = "\t",
    col.names = TRUE
)