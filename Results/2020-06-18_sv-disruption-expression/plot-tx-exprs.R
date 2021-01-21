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
                est_counts / eff_len / sf + offset,
                base = base
            )
        ),
        keyby = c("sample", "target_id", "transcript_name")
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
so <- readRDS(file.path("sleuth", "test_201.sleuth-object.rds"))

mut_samples <- "PCa53687"

nonmut_samples <- setdiff(
    c(
        "PCa13266", "PCa13848", "PCa14121", "PCa19121", "PCa3023", "PCa33173",
        "PCa40507", "PCa51852", "PCa53687", "PCa56413", "PCa57294", "PCa58215"
    ),
    mut_samples
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

ALT_GENES <- names(so)

braf_tx_ids <- gencode_tx[gene_name == "BRAF", ens_transcript]
tx_counts <- as.data.table(so$obs_norm)[target_id %in% braf_tx_ids]
scale_factors <- data.table(
    sample = names(so$est_counts_sf),
    sf = so$est_counts_sf
)
tx_counts <- merge(
    x = tx_counts,
    y = scale_factors,
    by = "sample"
)

# set order of transcripts by mean TPM of non-mutant samples
tx_order <- tx_counts[sample != mut_sample, mean(tpm), by = target_id][, target_id, keyby = V1]$target_id
tx_counts[, target_id := factor(target_id, ordered = TRUE, levels = rev(tx_order))]

# ==============================================================================
# Plots
# ==============================================================================
loginfo("Plotting data")

gg <- (
    ggplot()
    + geom_point(
        data = tx_counts[sample != mut_samples],
        mapping = aes(
            x = target_id,
            y = tpm
        ),
        fill = "#BDBDBD",
        colour = "black",
        position = position_jitter(height = 0, width = 0.3),
        shape = 21,
        size = 4,
        alpha = 0.5
    )
    + geom_point(
        data = tx_counts[sample == mut_sample],
        mapping = aes(
            x = target_id,
            y = tpm
        ),
        fill = "#ff6347",
        colour = "black",
        shape = 21,
        size = 4
    )
    + geom_boxplot(
        data = tx_counts[sample %in% nonmut_samples],
        mapping = aes(
            x = target_id,
            y = tpm
        ),
        outlier.shape = NA,
        alpha = 0.7
    )
    + labs(
        x = "BRAF Transcript",
        y = "TPM"
    )
    + theme_minimal()
    + theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    )
)
savefig(gg, "Plots/tx-exprs.BRAF", height = 10, width = 6)
