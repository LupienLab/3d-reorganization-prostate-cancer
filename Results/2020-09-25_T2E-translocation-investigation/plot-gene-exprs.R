# ==============================================================================
# Meta
# ==============================================================================
# Plot gene expression
# --------------------------------------
# Description: Plot the gene expression changes for the T2E translocation insertion site
# Author: James Hawley


# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("logging"))

loginfo("Loading packages")
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
suppressMessages(library("sleuth"))

ALT_STATUS <- c("altered", "notaltered")

# ==============================================================================
# Functions
# ==============================================================================
#' Aggregate estimated transcript counts into an estimated gene count
# 
#' Description
#'
#' @param x Estimated transcript counts
#' @param x Estimated transcript lengths
#' @param sf Sample scaling factor
#' @param offset Logarithm offset
agg_tx <- function(x, lens, sf, offset = 0.5, base = exp(1)) {
    med_length <- median(lens)
    # if sf is a vector, only take the first element
    sf <- sf[1]
    return(log(med_length / sf * sum(x / lens) + offset, base = base))
}

#' Save figures in multiple formats
#'
#' @param gg ggplot object
#' @param prefix Prefix for output file
#' @param ext Output extensions
#' @param dpi DPI resolution
savefig <- function(gg, prefix, ext = c("png", "pdf"), width = 20, height = 12, dpi = 400) {
    for (e in ext) {
        ggsave(
            paste(prefix, e, sep = "."),
            gg,
            height = height,
            width = width,
            units = "cm",
            dpi = dpi
        )
    }
}

# ==============================================================================
# Data
# ==============================================================================
loginfo("Loading data")

# load Kallisto/Sleuth object
so <- list(
    "altered" = readRDS(file.path("..", "2020-06-18_sv-disruption-expression", "sleuth", "test_44.sleuth-object.rds")),
    "notaltered" = readRDS(file.path("..", "2020-06-18_sv-disruption-expression", "sleuth", "test_197.sleuth-object.rds"))
)

mut_samples <- list(
    "altered" = "PCa13848",
    "notaltered" = "PCa53687"
)

# load GENCODE annotations
gencode_genes <- fread(
    file.path("..", "..", "Data", "External", "GENCODE", "gencode.v33.all-genes.bed"),
    sep = "\t",
    header = FALSE,
    col.names = c("chr", "start", "end", "strand", "ens_gene", "gene_name")
)
gencode_tx <- fread(
    file.path("..", "..", "Data", "External", "GENCODE", "gencode.v33.all-transcripts.bed"),
    sep = "\t",
    header = FALSE,
    col.names = c("chr", "start", "end", "strand", "ens_gene", "gene_name", "ens_transcript", "transcript_name")
)

# ==============================================================================
# Analysis
# ==============================================================================
loginfo("Performing calculations")

# the genes within the area of interest
test_genes <- list(
    "altered" = c("INSM2", "NFKBIA", "RALGAPA1", "BRMS1L"),
    "notaltered" = c(
        "SLC26A4",
        "COG5",
        "BCAP29",
        "SLC26A3",
        "DUS4L",
        "CBLL1",
        "WBP1LP2"
    )
)

# extract annotations for these genes
these_genes <- lapply(
    ALT_STATUS,
    function(alt_status) {
        gencode_genes[gene_name %in% test_genes[[alt_status]]]
    }
)
names(these_genes) <- ALT_STATUS
these_tx <- lapply(
    ALT_STATUS,
    function(alt_status) {
        gencode_tx[gene_name %in% test_genes[[alt_status]]]
    }
)
names(these_tx) <- ALT_STATUS

# get kallisto/sleuth normalized counts
so_counts <- lapply(
    ALT_STATUS,
    function(alt_status) {
        as.data.table(so[[alt_status]]$obs_norm)
    }
)
names(so_counts) <- ALT_STATUS
so_counts <- lapply(
    ALT_STATUS,
    function(alt_status) {
        so_counts[[alt_status]][target_id %in% these_tx[[alt_status]]$ens_transcript]
    }
)
names(so_counts) <- ALT_STATUS
so_counts <- lapply(
    ALT_STATUS,
    function(alt_status) {
        merge(
            x = these_tx[[alt_status]],
            y = so_counts[[alt_status]],
            by.x = "ens_transcript",
            by.y = "target_id"
        )
    }
)
names(so_counts) <- ALT_STATUS

# get scaling factors for each sample and merge into the table
sf <- lapply(
    ALT_STATUS,
    function(alt_status) {
        data.table(
            "SampleID" = names(so[[alt_status]]$est_counts_sf),
            "Scale" = so[[alt_status]]$est_counts_sf
        )
    }
)
names(sf) <- ALT_STATUS

so_counts <- lapply(
    ALT_STATUS,
    function(alt_status) {
        merge(
            x = sf[[alt_status]],
            y = so_counts[[alt_status]],
            by.x = "SampleID",
            by.y = "sample"
        )
    }
)
names(so_counts) <- ALT_STATUS

# calculate aggregated gene count
so_gene_counts <- lapply(
    ALT_STATUS,
    function(alt_status) {
        so_counts[[alt_status]][,
            lapply(.SD, function(row) agg_tx(est_counts, eff_len, Scale, 2)),
            keyby = c("SampleID", "gene_name", "ens_gene")
        ][, .SD, .SDcols = c("SampleID", "gene_name", "ens_gene", "est_counts")]
    }
)
names(so_gene_counts) <- ALT_STATUS

# ==============================================================================
# Plots
# ==============================================================================
loginfo("Plotting data")

gg <- list()

for (alt_status in c("altered", "notaltered")) {
    gg[[alt_status]] <- (
        ggplot()
        + geom_boxplot(
            data = so_gene_counts[[alt_status]][SampleID != mut_samples[[alt_status]]],
            mapping = aes(x = gene_name, y = est_counts),
            outlier.shape = NA
        )
        + geom_point(
            data = so_gene_counts[[alt_status]][SampleID != mut_samples[[alt_status]]],
            mapping = aes(x = gene_name, y = est_counts),
            colour = "#BDBDBD",
            position = position_jitter(height = 0, width = 0.3)
        )
        + geom_point(
            data = so_gene_counts[[alt_status]][SampleID == mut_samples[[alt_status]]],
            mapping = aes(x = gene_name, y = est_counts),
            colour = "#ff6347"
        )
        + labs(x = "Gene", y = expression("Estimated mRNA count (" * log[2] * ")"))
        + theme_minimal()
        + theme(
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
        )
    )
    savefig(gg[[alt_status]], paste0("gene-exprs.", alt_status, "-TAD"), height = 10, width = 6)
}
