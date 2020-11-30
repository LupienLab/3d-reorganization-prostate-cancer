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
    "ERG" = readRDS(file.path("sleuth", "test_46.sleuth-object.rds")),
    "TMPRSS2" = readRDS(file.path("sleuth", "test_45.sleuth-object.rds")),
    "ZNF516" = readRDS(file.path("sleuth", "test_42.sleuth-object.rds")),
    "PMEPA1" = readRDS(file.path("sleuth", "test_43.sleuth-object.rds")),
    "BRAF" = readRDS(file.path("sleuth", "test_201.sleuth-object.rds")),
    "DENND2A" = readRDS(file.path("sleuth", "test_201.sleuth-object.rds")),
)

mut_samples <- list(
    "ERG" = c("PCa13848", "PCa19121", "PCa3023", "PCa40507", "PCa51852", "PCa56413"),
    "TMPRSS2" = c("PCa13848", "PCa19121", "PCa3023", "PCa40507", "PCa51852", "PCa56413"),
    "ZNF516" = c("PCa13848"),
    "PMEPA1" = c("PCa13848"),
    "BRAF" = c("PCa53687"),
    "DENND2A" = c("PCa53687")
)
nonmut_samples <- lapply(
    names(mut_samples),
    function(alt_gene) {
        return(setdiff(
            c(
                "PCa13266", "PCa13848", "PCa14121", "PCa19121", "PCa3023", "PCa33173",
                "PCa40507", "PCa51852", "PCa53687", "PCa56413", "PCa57294", "PCa58215"
            ),
            mut_samples[[alt_gene]]
        ))
    }
)
names(nonmut_samples) <- names(mut_samples)

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

ALT_GENES <- names(so)

# extract annotations for these genes
these_genes <- gencode_genes[gene_name %in% ALT_GENES]
these_tx <- lapply(
    ALT_GENES,
    function(alt_gene) {
        gencode_tx[gene_name == alt_gene]
    }
)
names(these_tx) <- ALT_GENES

# get kallisto/sleuth normalized counts
so_counts <- lapply(
    ALT_GENES,
    function(alt_gene) {
        as.data.table(so[[alt_gene]]$obs_norm)
    }
)
names(so_counts) <- ALT_GENES

so_counts <- lapply(
    ALT_GENES,
    function(alt_gene) {
        so_counts[[alt_gene]][target_id %in% these_tx[[alt_gene]]$ens_transcript]
    }
)
names(so_counts) <- ALT_GENES
so_counts <- lapply(
    ALT_GENES,
    function(alt_gene) {
        merge(
            x = these_tx[[alt_gene]],
            y = so_counts[[alt_gene]],
            by.x = "ens_transcript",
            by.y = "target_id"
        )
    }
)
names(so_counts) <- ALT_GENES

# get scaling factors for each sample and merge into the table
sf <- lapply(
    ALT_GENES,
    function(alt_gene) {
        data.table(
            "SampleID" = names(so[[alt_gene]]$est_counts_sf),
            "Scale" = so[[alt_gene]]$est_counts_sf
        )
    }
)
names(sf) <- ALT_GENES

so_counts <- lapply(
    ALT_GENES,
    function(alt_gene) {
        merge(
            x = sf[[alt_gene]],
            y = so_counts[[alt_gene]],
            by.x = "SampleID",
            by.y = "sample"
        )
    }
)
names(so_counts) <- ALT_GENES

# calculate aggregated gene count
so_gene_counts <- lapply(
    ALT_GENES,
    function(alt_gene) {
        so_counts[[alt_gene]][,
            lapply(.SD, function(row) agg_tx(est_counts, eff_len, Scale, 2)),
            keyby = c("SampleID", "gene_name", "ens_gene")
        ][, .SD, .SDcols = c("SampleID", "gene_name", "ens_gene", "est_counts")]
    }
)
names(so_gene_counts) <- ALT_GENES

# ==============================================================================
# Plots
# ==============================================================================
loginfo("Plotting data")

gene_pairs <- list(
    c("ERG", "TMPRSS2")
    #c("ZNF516", "PMEPA1"),
    #c("BRAF", "DENND2A")
)

for (gp in gene_pairs) {
    dt <- rbindlist(so_gene_counts[gp])
    # combined mut samples
    ms <- unique(unlist(mut_samples[gp]))
    nms <- unique(unlist(nonmut_samples[gp]))
    gg <- (
        ggplot()
        + geom_boxplot(
            data = dt[SampleID %in% nms],
            mapping = aes(x = gene_name, y = est_counts),
            outlier.shape = NA
        )
        + geom_point(
            data = dt[SampleID %in% nms],
            mapping = aes(x = gene_name, y = est_counts),
            colour = "#BDBDBD",
            position = position_jitter(height = 0, width = 0.3)
        )
        + geom_point(
            data = dt[SampleID %in% ms],
            mapping = aes(x = gene_name, y = est_counts),
            colour = "#ff6347",
            position = position_jitter(height = 0, width = 0.3)
        )
        + labs(x = "Gene", y = expression("Estimated mRNA count (" * log[2] * ")"))
        + theme_minimal()
        + theme(
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
        )
    )
    savefig(gg, paste0("Plots/gene-exprs.", paste(gp, collapse = "-")), height = 10, width = 6)
}
