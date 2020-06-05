# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("Matrix"))
suppressMessages(library("data.table"))

if (!interactive()) {
    PARSER <- argparse::ArgumentParser(
        description = "Convert sparse COO matrix to wide dense matrix"
    )
    PARSER$add_argument(
        "mtx",
        type = "character",
        help = "Input matrix file"
    )
    PARSER$add_argument(
        "out",
        type = "character",
        help = "Output matrix file"
    )
    ARGS <- PARSER$parse_args()
} else {
    ARGS <- list(
        "mtx" = "TMP/PCa51852.chr7_137700000_139700000.mtx",
        "out" = "TMP/PCa51852.chr7_137700000_139700000.wide.mtx"
    )
}

# ==============================================================================
# Data
# ==============================================================================
# load matrix dump
mat <- fread(ARGS$mtx, sep = "\t", header = TRUE, col.names = c("i", "j", "v"))
# adjust indices to keep matrix small
id_offset <- mat[, min(i)]
mat <- mat[, .(i = i - id_offset + 1, j = j - id_offset + 1, v = v)]

# get dimension of matrix
n <- mat[, max(i, j)]

# convert to sparse matrix
mtx <- sparseMatrix(i = mat$i, j = mat$j, x = mat$v, symmetric = TRUE, dims = c(n, n))
mtx_dense <- as.matrix(mtx)

# ==============================================================================
# Save data
# ==============================================================================
write.table(mtx_dense, file = ARGS$out, row.names = FALSE, col.names = FALSE, sep = "\t")
