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

RES_DIR <- file.path(
    "..", "..", "results", "2020-06-18_sv-disruption-expression"
)
PLOT_DIR <- file.path(RES_DIR, "Plots")

# ==============================================================================
# Functions
# ==============================================================================
get_important_sleuth_info <- function(path, alt_gene, offset = 0.5, base = exp(1)) {
    print(alt_gene)
    # extract annotations for these genes
    gene <- gencode_genes[gene_name == alt_gene]
    gene_tx <- gencode_tx[gene_name == alt_gene]

    # get kallisto/sleuth normalized counts
    sleuth_obj <- readRDS(path)
    norm_counts <- as.data.table(sleuth_obj$obs_norm)
    # keep only the counts related to this gene and its transcripts
    tx_counts <- norm_counts[target_id %in% gene_tx$ens_transcript]
    # merge annotation with count information
    tx_counts <- merge(
        x = gene_tx,
        y = tx_counts,
        by.x = "ens_transcript",
        by.y = "target_id"
    )
    # get scaling factors for each sample and merge into the table
    sf <- data.table(
        "sample" = names(sleuth_obj$est_counts_sf),
        "sf" = sleuth_obj$est_counts_sf
    )
    tx_counts <- merge(
        x = tx_counts,
        y = sf,
        by = "sample"
    )
    gene_counts <- tx_counts[,
        .(
            est_counts = log(
                median(eff_len) / sf * sum(est_counts / eff_len) + offset,
                base = base
            )
        ),
        keyby = c("sample", "ens_gene", "gene_name")
    ]
    return(list(
        "tx" = tx_counts,
        "gene" = gene_counts
    ))
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
    "ERG" = file.path(RES_DIR, "sleuth", "test_46.sleuth-object.rds"),
    "TMPRSS2" = file.path(RES_DIR, "sleuth", "test_45.sleuth-object.rds"),
    "ZNF516" = file.path(RES_DIR, "sleuth", "test_42.sleuth-object.rds"),
    "PMEPA1" = file.path(RES_DIR, "sleuth", "test_43.sleuth-object.rds")
)

mut_samples <- list(
    "ERG" = c("PCa13848", "PCa19121", "PCa3023", "PCa40507", "PCa51852", "PCa56413"),
    "TMPRSS2" = c("PCa13848", "PCa19121", "PCa3023", "PCa40507", "PCa51852", "PCa56413"),
    "ZNF516" = c("PCa13848"),
    "PMEPA1" = c("PCa13848")
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
    file.path("..", "..", "data", "External", "GENCODE", "gencode.v33.all-genes.bed"),
    sep = "\t",
    header = FALSE,
    col.names = c("chr", "start", "end", "strand", "ens_gene", "gene_name")
)
gencode_tx <- fread(
    file.path("..", "..", "data", "External", "GENCODE", "gencode.v33.all-transcripts.bed"),
    sep = "\t",
    header = FALSE,
    col.names = c("chr", "start", "end", "strand", "ens_gene", "gene_name", "ens_transcript", "transcript_name")
)

# ==============================================================================
# Analysis
# ==============================================================================
loginfo("Performing calculations")

ALT_GENES <- names(so)

# calculate aggregated gene count
so_counts <- lapply(
    ALT_GENES,
    function(alt_gene) {
        get_important_sleuth_info(so[[alt_gene]], alt_gene)
    }
)
names(so_counts) <- ALT_GENES

# ==============================================================================
# Plots
# ==============================================================================
loginfo("Plotting data")

gene_pairs <- list(
    c("ERG", "TMPRSS2"),
    c("ZNF516", "PMEPA1")
)

for (gp in gene_pairs) {
    dt <- unique(rbindlist(lapply(gp, function(alt_gene) so_counts[[alt_gene]]$gene)))
    # combined mut samples
    ms <- unique(unlist(mut_samples[gp]))
    nms <- unique(unlist(nonmut_samples[gp]))
    gg <- (
        ggplot()
        + geom_boxplot(
            data = dt[sample %in% nms],
            mapping = aes(x = gene_name, y = est_counts),
            outlier.shape = NA
        )
        + geom_point(
            data = dt[sample %in% nms],
            mapping = aes(x = gene_name, y = est_counts),
            colour = "#BDBDBD",
            position = position_jitter(height = 0, width = 0.3)
        )
        + geom_point(
            data = dt[sample %in% ms],
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
    savefig(
        gg,
        file.path(
            PLOT_DIR, 
            paste0("gene-exprs.", paste(gp, collapse = "-"))
        ),
        height = 10,
        width = 6
    )
}
