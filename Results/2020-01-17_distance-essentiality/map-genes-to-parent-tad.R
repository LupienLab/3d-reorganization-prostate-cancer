# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("argparse"))

if (!interactive()) {
    PARSER <- argparse::ArgumentParser(
        description = "Calculate the distribution of genes across the length of TADs"
    )
    PARSER$add_argument(
        "TADs",
        type = "character",
        help = "Aggregated TAD calls"
    )
    PARSER$add_argument(
        "nearest",
        type = "character",
        help = "Boundaries nearest to each gene's TSS"
    )
    PARSER$add_argument(
        "-r", "--resolution",
        type = "integer",
        help = "Contact matrix resolution in bp",
        default = 40000
    )
    PARSER$add_argument(
        "-p", "--prefix",
        type = "character",
        help = "Prefix for output files",
        default = "output"
    )
    ARGS <- PARSER$parse_args()
} else {
    ARGS = list(
        TADs = "../2020-01-15_TAD-aggregation/resolved-TADs/PCa3023.40000bp.aggregated-domains.sorted.bedGraph",
        nearest = "Closest/PCa3023.closest-boundaries.bed",
        resolution = 40000,
        prefix = "Closest/PCa3023"
    )
}

# ==============================================================================
# Data
# ==============================================================================
# load TADs
tads = fread(
    ARGS$TADs,
    sep = "\t",
    header = FALSE,
    col.names = c("chr", "start", "end", "order_lower", "order_upper", "w")
)

# load genes and their nearest boundaries
# `w` column is a list of integers, separated by "|"
genes = fread(
    ARGS$nearest,
    sep = "\t",
    sep2 = "|", # not implemented in data.table yet
    header = FALSE,
    col.names = c(
        "chr", "start_TSS", "end_TSS", "name", "score", "strand", "Ensembl_ID",
        "chr_bound", "start_bound", "end_bound", "order", "w", "distance"
    )
)


# ==============================================================================
# Analysis
# ==============================================================================
# figure out which TAD each gene belongs to
## I know for loops aren't the most efficient in R, but there are ~20K genes, so this is fine enough
genes_copy = copy(genes)
genes = genes_copy
for (i in 1:genes[, .N]) {
    if (i %% 1000 == 0) cat(".")
    # get the smallest TAD that encapsulates this gene
    parent_tads = tads[chr == genes[i, chr] & start <= genes[i, start_TSS] & end >= genes[i, end_TSS]]
    # if this gene straddles domains across all length scales, remove it
    if (parent_tads[, .N] == 0) {
        genes[i, start_TAD := -1]
        genes[i, end_TAD := -1]
        next
    }
    # select the parent TAD with a boundary that is closest to the gene
    closest_parent_tad = parent_tads[start == genes[i, start_bound] | end == genes[i, start_bound], .SD]
    # if none of the parents match the closest boundary, check +/- a bin to resolve this
    # I haven't run into this situation, but I'm going to leave this in just in case
    if (closest_parent_tad[, .N] == 0 ) {
        # if (closest_parent_tad$start != genes[i, start_TAD] & closest_parent_tad$end != genes[i, start_TAD] & abs(genes[i, distance]) < ARGS$resolution) {
        # } else {
        # }
        cat("Gene\n")
        print(genes[i, .(chr, start, end, name, strand, start_bound, end_bound, order)])
        cat("\nAll Parents\n")
        print(parent_tads)
        cat("\nClosest Parent\n")
        print(closest_parent_tad)
        dn_near_tads = tads[chr == genes[i, chr] & end == genes[i, start_bound] & order_upper == genes[i, order]]
        up_near_tads = tads[chr == genes[i, chr] & start == genes[i, start_bound] & order_lower == genes[i, order]]
        cat("\nNearby\n")
        print(dn_near_tads)
        print(up_near_tads)
        cat("\n\n")
        Sys.sleep(2)
    # if there's no singular TAD with the closest boundary (possible with calls at multiple windows)
    } else if (closest_parent_tad[, .N] > 1) {
        # pick the smallest in size
        # not using which.min here, since multiple w's are possible and this returns the first
        closest_parent_tad = closest_parent_tad[(end - start) == min(end - start), .SD]
        # pick the TAD with the smallest size, then smallest window size if necessary
        closest_parent_tad = closest_parent_tad[which.min(w), .SD]
    }
    if (closest_parent_tad[, .N] > 1) print("help") # I haven't run into this situation, but I'm going to leave this in just in case
    # assign a parent TAD to each gene, along with all of that TADs information
    genes[i, start_TAD := closest_parent_tad$start]
    genes[i, end_TAD := closest_parent_tad$end]
    genes[i, w_TAD := closest_parent_tad$w]
    genes[i, order_lower_TAD := closest_parent_tad$order_lower]
    genes[i, order_upper_TAD := closest_parent_tad$order_upper]
}
# record what percentage of the distance through the TAD the beginning of each gene is
# ignore genes that straddle boundaries at all window sizes
genes = genes[start_TAD != -1]
genes[, Fraction := ifelse(
    strand == "+",
    0.5 - abs(start_TSS - (start_TAD + end_TAD) / 2) / (end_TAD - start_TAD),
    0.5 - abs(end_TSS - (start_TAD + end_TAD) / 2) / (end_TAD - start_TAD)
)]

# ==============================================================================
# Save data
# ==============================================================================
fwrite(
    # only write relevant information
    genes[, .SD, .SDcols = c(
        "chr", "start_TSS", "end_TSS", "name", "score", "strand", "Ensembl_ID",
        "start_TAD", "end_TAD", "w_TAD", "order_lower_TAD", "order_upper_TAD", "Fraction"
    )],
    paste0(ARGS$prefix, ".distance-dependency.tsv"),
    sep = "\t",
    col.names = TRUE
)
