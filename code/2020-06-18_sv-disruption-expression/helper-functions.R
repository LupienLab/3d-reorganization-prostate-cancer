#' Get estimated mRNA counts for genes and transcripts from a specified file
# '
#' @param path Path to sleuth object
#' @param alt_gene Gene name
#' @param offset Default offset value inside the logarithm for stabilization
#' @param base Default logarithm base
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
        x = tx_counts,
        y = gene_tx,
        by.x = "target_id",
        by.y = "ens_transcript"
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
    tx_counts_norm <- tx_counts[,
        .(
            est_counts = log(
                est_counts / eff_len / sf + offset,
                base = base
            )
        ),
        keyby = c("sample", "target_id", "transcript_name")
    ]
    gene_counts_norm <- unique(tx_counts[,
        .(
            est_counts = log(
                median(eff_len) / sf * sum(est_counts / eff_len) + offset,
                base = base
            )
        ),
        keyby = c("sample", "ens_gene", "gene_name")
    ])

    # get differential analysis results
    # (all objects processed in this folder have this structure)
    so_genes <- as.data.table(sleuth_results(
        sleuth_obj,
        "conditionMutated",
        "wt",
        show_all = FALSE,
        pval_aggregate = TRUE
    ))
    # rows are duplicated with transcript IDs, but everything else is the same => unique only keeps relevant info
    so_genes <- unique(so_genes[,
        .(
            target_id,
            gene_name,
            num_aggregated_transcripts,
            sum_mean_obs_counts,
            pval,
            qval
        )
    ])

    so_transcripts <- as.data.table(sleuth_results(
        sleuth_obj,
        "conditionMutated",
        "wt",
        show_all = FALSE,
        pval_aggregate = FALSE
    ))
    # NAs in start/end positions force the columns to be "character" instead of "integer"
    # this happens because of the target mapping to GENCODE
    # I've only included canonical chromosomes, not extra scaffolds, but the original index that kallisto quantifies against
    # will contain these scaffolds, resulting in an NA target mapping from sleuth
    so_transcripts <- so_transcripts[complete.cases(so_transcripts)]
    so_transcripts[, `:=`(start = as.integer(start), end = as.integer(end))]

    return(list(
        "tx_counts" = tx_counts_norm,
        "gene_counts" = gene_counts_norm,
        "tx_de" = so_transcripts[gene_name == alt_gene],
        "gene_de" = so_genes[gene_name == alt_gene]
    ))
}
