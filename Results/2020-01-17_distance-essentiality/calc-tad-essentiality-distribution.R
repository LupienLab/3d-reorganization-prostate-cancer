# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
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

# calculate mean essentiality by TAD
# because some genes span TADs, they intersections needs to be double counted
# thus double counting genes, but there are only 1331 cases here, not enough
# to change the distribution overall of ~ 16K genes with essentiality scores
ess_overlaps = ess_coords[hits$queryHits, .SD]
ess_overlaps[, TAD := hits$subjectHits]

ess_long = melt(
    ess_overlaps,
    id.vars = c("Gene", "TAD", "EnsemblID"),
    measure.vars = c("22Rv1", "DU145", "LNCaP", "MDA PCa 2b", "NCI-H660", "PC3", "VCaP"),
    variable.name = "CellLine",
    value.name = "Essentiality"
)

tad_essentiality = ess_long[, .(mean(Essentiality, na.rm = TRUE), sd(Essentiality, na.rm = TRUE)), by = c("CellLine", "TAD")]

# ==============================================================================
# Plots
# ==============================================================================
gg = (
    ggplot(data = tad_essentiality)
    + geom_density(aes(x = V1))
    + labs(x = "Mean TAD Essentiality", y = "Density")
    + scale_x_continuous(
        limits = c(-2, 2)
    )
    + facet_wrap(~ CellLine)
    + theme_minimal()
)
ggsave(
    "Plots/TAD-mean-essentiality-density.png",
    height = 12,
    width = 20,
    units = "cm"
)

gg = (
    ggplot(data = tad_essentiality)
    + geom_density(aes(x = V2))
    + labs(x = "Mean TAD Essentiality", y = "Density")
    + scale_x_continuous(
        limits = c(-2, 2)
    )
    + facet_wrap(~ CellLine)
    + theme_minimal()
)
ggsave(
    "Plots/TAD-sd-essentiality-density.png",
    height = 12,
    width = 20,
    units = "cm"
)
