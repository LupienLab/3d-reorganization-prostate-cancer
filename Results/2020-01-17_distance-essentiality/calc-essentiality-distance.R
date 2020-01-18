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
        "-p", "--prefix",
        type = "character",
        help = "Prefix for output files",
        default = "output"
    )
    ARGS <- PARSER$parse_args()
} else {
    ARGS = list(
        TADs = "../2020-01-15_TAD-aggregation/resolved-TADs/PCa3023.40000bp.aggregated-domains.sorted.bedGraph",
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
    file.path("..", "..", "Data", "External", "GENCODE", "gencode.v33.genes.sorted.bed"),
    sep = "\t",
    sep2 = "|", # not implemented in data.table yet
    header = FALSE,
    drop = 5, # drop score column
    col.names = c(
        "chr", "start", "end", "name", "strand", "Ensembl_ID"
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
    # a local copy of the row (changes to curr_gene != changes to genes[i])
    curr_gene = genes[i]
    # get the smallest TAD that encapsulates this gene
    parent_tads = tads[chr == curr_gene$chr & start <= curr_gene$start & end >= curr_gene$end]
    # if this gene straddles domains across all length scales, remove it
    if (parent_tads[, .N] == 0) {
        genes[i, start_tad := -1]
        genes[i, end_tad := -1]
    } else {
        smallest_parent_tad = parent_tads[which.min(w), .SD]
        genes[i, start_tad := smallest_parent_tad$start]
        genes[i, end_tad := smallest_parent_tad$end]
    }
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
