# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("argparse"))

if (!interactive()) {
    PARSER <- argparse::ArgumentParser(
        description = "Calculate the distribution of essentiality scores of genes across the length of TADs"
    )
    PARSER$add_argument(
        "TADs",
        type = "character",
        help = "Aggregated TAD calls"
    )
    PARSER$add_argument(
        "bounds",
        type = "character",
        help = "Boundaries nearest to each gene's TSS"
    )
    PARSER$add_argument(
        "-r",, "--resolution",
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
        prefix = "essential-distance/PCa3023"
    )
}

# ==============================================================================
# Data
# ==============================================================================
# load DepMap essentiality data
depmap = fread(file.path("..", "..", "Data", "External", "DepMap", "depmap-rnai.tsv"))

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
        "chr", "start", "end", "name", "score", "strand", "Ensembl_ID",
        "chr_TAD", "start_TAD", "end_TAD", "order", "w", "distance"
    )
)

# map essentiality scores to genes
genes = merge(
    x = genes,
    y = depmap,
    by.x = "name",
    by.y = "Gene",
    all.x = TRUE,
    all.y = FALSE
)

# ==============================================================================
# Analysis
# ==============================================================================
# figure out which TAD each gene belongs to
genes_copy = copy(genes)
genes = genes_copy
for (i in 1:genes[, .N]) {
    if (i %% 1000 == 0) cat(".")
    # get the smallest TAD that encapsulates this gene
    parent_tads = tads[chr == genes[i, chr] & start <= genes[i, start] & end >= genes[i, end]]
    # if this gene straddles domains across all length scales, remove it
    if (parent_tads[, .N] == 0) {
        genes[i, start_tad := -1]
        genes[i, end_tad := -1]
        next
    }
    # select the parent TAD with a boundary that is closest to the gene
    closest_parent_tad = parent_tads[start == genes[i, start_TAD] | end == genes[i, start_TAD], .SD]
    # if none of the parents match the closest boundary, check +/- a bin to resolve this
    if (closest_parent_tad[, .N] == 0 ) {
        # if (closest_parent_tad$start != genes[i, start_TAD] & closest_parent_tad$end != genes[i, start_TAD] & abs(genes[i, distance]) < ARGS$resolution) {
        # } else {
        # }
        cat("Gene\n")
        print(genes[i, .(chr, start, end, name, strand, start_TAD, end_TAD, order)])
        cat("\nAll Parents\n")
        print(parent_tads)
        cat("\nClosest Parent\n")
        print(closest_parent_tad)
        dn_near_tads = tads[chr == genes[i, chr] & end == genes[i, start_TAD] & order_upper == genes[i, order]]
        up_near_tads = tads[chr == genes[i, chr] & start == genes[i, start_TAD] & order_lower == genes[i, order]]
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
        closest_parent_tad = closest_parent_tad[which.min(w), .SD]
    }
    if (closest_parent_tad[, .N] > 1) print("help")
    genes[i, start_tad := closest_parent_tad$start]
    genes[i, end_tad := closest_parent_tad$end]
}
# record what percentage of the distance through the TAD the beginning of each gene is
# ignore genes that straddle boundaries at all window sizes
genes = genes[start_tad != -1]
genes[, Fraction := ifelse(
    strand == "+",
    0.5 - abs(start - (start_tad + end_tad) / 2) / (end_tad - start_tad),
    0.5 - abs(end - (start_tad + end_tad) / 2) / (end_tad - start_tad)
)]

# ==============================================================================
# Save data
# ==============================================================================
fwrite(
    genes,
    paste0(ARGS$prefix, ".distance-dependency.tsv"),
    sep = "\t",
    col.names = TRUE
)
