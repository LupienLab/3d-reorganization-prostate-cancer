# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("sleuth"))
suppressMessages(library("data.table"))

KALLISTO_DIR <- file.path("..", "..", "Data", "Processed", "2020-06-17_PCa-RNA-seq")

# ==============================================================================
# Data
# ==============================================================================
# load sample meta
meta <- fread(
    "config.tsv",
    sep = "\t",
    header = TRUE
)
meta <- meta[(Include == "Yes") & (Type == "Malignant")]
SAMPLES <- meta$Sample_ID

# add kallisto paths
meta[, Path := file.path(KALLISTO_DIR, Sample_ID)]

# create experimental design matrix
design <- meta[, .(
    sample = Sample_ID,
    path = Path
)]

cmpmt_samples <- list(
    "benign_like" = c(
        "PCa13266",
        "PCa13848",
        "PCa3023",
        "PCa33173",
        "PCa40507",
        "PCa53687",
        "PCa56413",
        "PCa58215"
    ),
    "novel" = c(
        "PCa14121",
        "PCa19121",
        "PCa51852",
        "PCa57294"
    )
)
design[sample %in% cmpmt_samples[["benign_like"]], condition := "benign"]
design[sample %in% cmpmt_samples[["novel"]], condition := "novel"]

# load GENCODE annotations
transcripts <- fread(
    file.path("..", "..", "Data", "External", "GENCODE", "gencode.v33.all-transcripts.bed"),
    sep = "\t",
    header = FALSE,
    col.names = c("chr", "start", "end", "strand", "ens_gene", "gene_name", "target_id", "transcript_name")
)

# ==============================================================================
# Analysis
# ==============================================================================
# create sleuth object, aggregating transcripts to genes
so <- sleuth_prep(
    design,
    extra_bootstrap_summary = TRUE,
    target_mapping = transcripts,
    aggregation_column = "ens_gene",
    num_cores = 4
)

# fit full model
so <- sleuth_fit(so, ~condition, "full")

# perform differential analysis
so <- sleuth_wt(so, "conditionnovel")

# extract results, focusing on chr19
so_transcripts <- as.data.table(sleuth_results(
    so,
    "conditionnovel",
    "wt",
    show_all = FALSE,
    pval_aggregate = FALSE
))[chr == "chr19"]

# re-calculate q-value since we only care about chr19 in this case
so_transcripts[, qval := p.adjust(pval, method = "fdr")]

# NAs in start/end positions force the columns to be "character" instead of "integer"
# this happens because of the target mapping to GENCODE
# I've only included canonical chromosomes, not extra scaffolds, but the original index that kallisto quantifies against
# will contain these scaffolds, resulting in an NA target mapping from sleuth
so_transcripts[, `:=`(start = as.integer(start), end = as.integer(end))]

so_genes <- as.data.table(sleuth_results(
    so,
    "conditionnovel",
    "wt",
    show_all = FALSE,
    pval_aggregate = TRUE
))[chr == "chr19"]

# rows are duplicated with transcript IDs, but everything else is the same => unique only keeps relevant info
so_genes <- unique(so_genes[, .(target_id, gene_name, num_aggregated_transcripts, sum_mean_obs_counts, pval, qval)])
# re-calculate q-value since we only care about chr19 in this case
so_genes[, qval := p.adjust(pval, method = "fdr")]


# ==============================================================================
# Save data
# ==============================================================================
# save annotated table
fwrite(so_genes, "dge.chr19-genes.tsv", sep = "\t", col.names = TRUE)
fwrite(so_transcripts, "dge.chr19-transcripts.tsv", sep = "\t", col.names = TRUE)
saveRDS(so, "dge.sleuth-object.rds")
