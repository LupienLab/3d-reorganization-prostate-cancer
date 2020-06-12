# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("CREAM"))
suppressMessages(library("regioneR"))
suppressMessages(library("argparse"))

if (!interactive()) {
    PARSER <- argparse::ArgumentParser(
        description = "Identify Clusters Of Regulatory Elements using CREAM"
    )
    PARSER$add_argument(
        "peaks",
        type = "character",
        help = "Input peaks file"
    )
    PARSER$add_argument(
        "output",
        type = "character",
        help = "Output file"
    )
    ARGS <- PARSER$parse_args()
} else {
    ARGS <- list(
        "input" = file.path("Peaks", "Pca13266_peaks.filtered.narrowPeak"),
        "output" = file.path("COREs", "PCa13266.cores.bed")
    )
}

# ==============================================================================
# Analysis
# ==============================================================================
# identify COREs
cores <- as.data.table(CREAM(in_path = ARGS$peaks))
colnames(cores) <- c("chr", "start", "end")

# convert to GRanges object
# BED files are 0-indexed but GRanges are 1-indexed
# I can safely ignore this offset here since I'm just sorting and writing, not calculating anything
cores_sorted <- as.data.table(sort(toGRanges(cores, genome = "hg38")))
colnames(cores_sorted)[1] <- "chr"

# ==============================================================================
# Save data
# ==============================================================================
fwrite(cores[, .SD, .SDcols = c("chr", "start", "end")], ARGS$output, sep = "\t", col.names = FALSE)
