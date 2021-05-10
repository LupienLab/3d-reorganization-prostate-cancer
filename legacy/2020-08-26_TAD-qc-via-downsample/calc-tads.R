# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("Matrix"))
suppressMessages(library("bit64"))
suppressMessages(library("argparse"))
source("topdom.R")

if (!interactive()) {
    PARSER <- argparse::ArgumentParser(
        description = "Call TADs on a contact matrix using TopDom"
    )
    PARSER$add_argument(
        "data",
        type = "character",
        help = "Annotated contact matrix file in (i, j, v) format"
    )
    PARSER$add_argument(
        "-w", "--window",
        type = "integer",
        help = "Window parameter for TopDom",
        default = 5
    )
    PARSER$add_argument(
        "-p", "--prefix",
        type = "character",
        help = "Prefix for output files",
        default = "contacts"
    )
    ARGS <- PARSER$parse_args()
} else {
    ARGS = list(
        data = "TMP/PCa13266.2000bp.chr2.annotated.mtx",
        prefix = "TADs/test",
        window = 3
    )
}

CHRS = paste0("chr", c(1:22, "X", "Y"))


# ==============================================================================
# Data
# ==============================================================================
# read cooler dump (only contains upper triangular elements)
vals = fread(ARGS$data, header = TRUE, sep = "\t")

# matrix IDs are from the genome-wide matrix, so to look at chromosome-specific
# bins, remove all the previous bins that aren't being measure
# i.e. shrink sparse matrix from N x N to M x M (M < N) without removing information
# about this particular chromosome, since all values in the first N - M pixels are 0
#   the (i, j, v) format is upper triangular, so j >= i for all rows
min_idx = vals[, min(bin1_id)]
max_idx = vals[, max(bin2_id)]
vals[, bin1_id := bin1_id - min_idx + 1]
vals[, bin2_id := bin2_id - min_idx + 1]

# dimension of the matrix
M = max_idx - min_idx + 1

# `cooler dump` either gives bins for the entire region or the matrix/pixel values
# but not all of the bins that matrix spans.
# TopDom requires the bins and the matrix to be the same dimension
# i.e. if the contact matrix is M x M, the `bins` data.frame should be M x 3

# get matrix bin size and chromosome
res = vals[1, end1 - start1]
chrom = vals[1, chrom1]

# get starting position of most 5' bin
#   unique just in case there are other non-zero pixels in the same row
start_pos = vals[bin1_id == 1, unique(start1)]

# get ending position of most 3' bin
#   unique just in case there are other non-zero pixels in the same column
end_pos = vals[bin2_id == M, unique(end2)]

# make contiguous bins
breakpoints = seq(start_pos, by = res, length.out= M + 1)
# adjust bins for end position of chromosome, if needed
breakpoints[M + 1] = min(end_pos, breakpoints[M + 1])

bins = data.table(
    chr = chrom,
    from.coord = breakpoints[1:M],
    to.coord = breakpoints[2:(M + 1)]
)

# remove rows with NAs in `balanced`
# vals = vals[!is.na(balanced), .SD]

# create sparse matrix from values
mat = sparseMatrix(
    i = vals[, bin1_id],
    j = vals[, bin2_id],
    x = vals[, balanced],
    symmetric = TRUE
)

# run TopDom to call TADs and boundaries
calls = TopDom(matrix.data = as.matrix(mat), bins = bins, window.size = ARGS$window)
# calls = TopDom(matrix.data = mat, bins = bins, window.size = ARGS$window)
# calls = TopDom(matrix.data = mat[1:50000, 1:50000], bins = bins[1:50000], window.size = ARGS$window, verbose = TRUE)

# save domains
bins_signal = as.data.table(calls$binSignal)
domains = as.data.table(calls$domain)
domains_bed = as.data.table(calls$bed)

fwrite(
    bins_signal,
    paste(ARGS$prefix, "bins-signal", "tsv", sep = "."),
    sep = "\t",
    col.names = TRUE
)

fwrite(
    domains,
    paste(ARGS$prefix, "domains", "tsv", sep = "."),
    sep = "\t",
    col.names = TRUE
)

fwrite(
    domains_bed,
    paste(ARGS$prefix, "domains", "bed", sep = "."),
    sep = "\t",
    col.names = FALSE
)

