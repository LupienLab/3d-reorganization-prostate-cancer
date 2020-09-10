# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("Matrix"))
suppressMessages(library("argparse"))



CHRS <- paste0("chr", c(1:22, "X"))

if (!interactive()) {
    PARSER <- argparse::ArgumentParser(
        description = "Merge contact matrices together"
    )
    PARSER$add_argument(
        "chrom",
        type = "character",
        help = "The chromosome to merge"
    )
    ARGS <- PARSER$parse_args()
} else {
    ARGS <- list(
        chrom = "chr22"
    )
}


# ==============================================================================
# Functions
# ==============================================================================
# simplify the process of reading in exported matrix files
read_contacts <- function(id, chrom) {
    # read zipped matrix file
    dt <- fread(
        cmd = paste0("zcat TMP/", id, ".300000000.res_40000bp.", chrom, ".noheader.mtx.gz"),
        sep = "\t"
    )
    # for some reason the last column is filled with NAs, but these shouldn't be here
    na_cols <- dt[, which(apply(.SD, 2, function(x) all(is.na(x))))]
    # remove these columns
    dt <- dt[, .SD, .SDcols = -na_cols]
    # convert to a matrix and remove column names
    mat <- as.matrix(dt)
    matnew <- Matrix(mat, sparse = TRUE)
    rm(mat)

    # return the matrix itself
    return(matnew)
}


# ==============================================================================
# Data
# ==============================================================================
cat("Loading Data\n")
# load metadata
metadata <- fread("config.tsv", sep = "\t", header = TRUE)
metadata <- metadata[Include == "Yes"]
SAMPLES <- list(
    "tumour" = metadata[Type == "Malignant", SampleID],
    "benign" = metadata[Type == "Benign", SampleID]
)


# ==============================================================================
# Analysis
# ==============================================================================

cat(ARGS$chrom, "\n")
for (i in c("tumour", "benign")) {
    cat("\t", i, "\n")
    matrices <- lapply(SAMPLES[[i]], read_contacts, ARGS$chrom)
    cat("\t\tSize of all sparse matrices in memory:", as.numeric(object.size(matrices)) / (1024 ^ 3), "GB\n")
    # calculate the mean contact matrix from all the tumour matrices
    mean_mat <- Reduce("+", matrices) / length(matrices)
    cat("\t\tSize of mean matrix in memory:", as.numeric(object.size(mean_mat)) / (1024 ^ 3), "GB\n")
    # convert to data.table for faster saving
    dt <- as.data.table(as.matrix(mean_mat))
    cat("\t\tSize of data.table:", as.numeric(object.size(dt)) / (1024 ^ 3), "GB\n")
    # save to output file
    fwrite(
        dt,
        paste0("TMP/", i, ".300000000.res_40000bp.", ARGS$chrom, ".mtx"),
        sep = "\t",
        col.names = FALSE
    )
    rm(dt, mean_mat, matrices)
}


