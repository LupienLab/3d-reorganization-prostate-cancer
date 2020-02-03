# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("argparse"))
suppressMessages(library("Matrix"))
suppressMessages(library("GenomicRanges"))

if (!interactive()) {
    PARSER <- argparse::ArgumentParser(
        description = "Calculate the distribution of essential genes across a set of TADs"
    )
    PARSER$add_argument(
        "tads",
        type = "character",
        help = "BED file of TADs"
    )
    ARGS <- PARSER$parse_args()
} else {
    ARGS = list(
        tads = "../2020-01-15_TAD-aggregation/resolved-TADs/separated-TADs/PCa13266.40000bp.w_30.domains.bed"
    )
}

# ==============================================================================
# Data
# ==============================================================================
# read TADs
tads = fread(ARGS$tads, sep = "\t", header = FALSE, col.names = c("chr", "start", "end", "lower_persistence", "upper_persistence"))

# read genes
genes = fread("../../Data/External/GENCODE/gencode.v33.genes.sorted.bed", sep = "\t", header = FALSE, col.names = c("chr", "start", "end", "name", "score", "strand", "EnsemblID"))

# read essentiality
essentiality = fread("../../Data/External/DepMap/depmap-rnai.tsv", sep = "\t", header = TRUE)


# ==============================================================================
# Analysis
# ==============================================================================
# merge essentiality and coordinate information
ess_coords = merge(
    x = essentiality,
    y = genes,
    by.x = "Gene",
    by.y = "name"
)

# create sparse matrix where each row is a TAD and each column is a percentage of that TAD
# i.e. an N x (1 / p) matrix where p \in [0, 1]
p = 0.1
mat = Matrix(0, nrow = tads[, .N], ncol = 1 / p, sparse = TRUE)

# find which genes map to which TADs
genes_gr = GRanges(
    seqnames = ess_coords$chr,
    ranges = IRanges(
        start = ess_coords$start + 1,
        end = ess_coords$end
    ),
    strand = ess_coords$strand,
    EnsemblID = ess_coords$EnsemblID,
    HUGO = ess_coords$Gene
)

tads_gr = GRanges(
    seqnames = tads$chr,
    ranges = IRanges(
        start = tads$start + 1,
        end = tads$end
    )
)

hits = as.data.table(findOverlaps(genes_gr, tads_gr))

# place genes into one of the 1/p bins
for (i in 1:tads[, .N]) {
    t = tads[i]
    t_gr = tads_gr[i]
    genes_in_tad = ess_coords[hits[subjectHits == i, queryHits], .SD]
    t_split = 
}

# ==============================================================================
# Plots
# ==============================================================================
