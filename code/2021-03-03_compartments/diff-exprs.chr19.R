# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("sleuth"))
suppressMessages(library("data.table"))

KALLISTO_DIR <- file.path("..", "..", "data", "Processed", "2020-06-17_PCa-RNA-seq")
RES_DIR <- file.path("..", "..", "results", "2021-03-03_compartments")
DGE_DIR <- file.path(RES_DIR, "DGE")

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
    file.path("..", "..", "data", "External", "GENCODE", "gencode.v33.all-transcripts.bed"),
    sep = "\t",
    header = FALSE,
    col.names = c("chr", "start", "end", "strand", "ens_gene", "gene_name", "target_id", "transcript_name")
)
genes <- fread(
    file.path("..", "..", "data", "External", "GENCODE", "gencode.v33.all-genes.bed"),
    sep = "\t",
    header = FALSE,
    col.names = c("chr", "start", "end", "strand", "ens_gene", "gene_name")
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

# extract results
so_transcripts <- as.data.table(sleuth_results(
    so,
    "conditionnovel",
    "wt",
    show_all = FALSE,
    pval_aggregate = FALSE
))

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
))

# rows are duplicated with transcript IDs, but everything else is the same => unique only keeps relevant info
so_genes <- unique(so_genes[, .(target_id, gene_name, num_aggregated_transcripts, sum_mean_obs_counts, pval, qval)])

# merge gene annotations
so_genes <- merge(
    x = so_genes,
    y = genes,
    by.x = c("target_id", "gene_name"),
    by.y = c("ens_gene", "gene_name"),
    all.x = TRUE
)

# simplify column order
so_genes <- so_genes[,
    .SD,
    keyby = "qval"
][,
    .SD,
    .SDcols = c(
        "chr", "start", "end", "strand",
        "target_id", "gene_name",
        "pval", "qval",
        "num_aggregated_transcripts", "sum_mean_obs_counts"
    )
]

# focus only on chr19 and chrY, where recurrent changes to compartmentalization are found
so_transcripts_cmpmt <- so_transcripts[(chr == "chr19") & (start > 20000000)]
so_genes_cmpmt <- so_genes[(chr == "chr19") & (start > 20000000)]

# re-calculate q-value since we only care about chr19 in this case
so_transcripts_cmpmt[, qval := p.adjust(pval, method = "fdr")]
so_genes_cmpmt[, qval := p.adjust(pval, method = "fdr")]


# ==============================================================================
# Save data
# ==============================================================================
# save annotated tables
fwrite(
    so_transcripts,
    file.path(DGE_DIR, "dge.chr19-diffs.all-transcripts.tsv"),
    sep = "\t",
    col.names = TRUE
)

fwrite(
    so_genes,
    file.path(DGE_DIR, "dge.chr19-diffs.all-genes.tsv"),
    sep = "\t",
    col.names = TRUE
)

fwrite(
    so_genes_cmpmt,
    file.path(DGE_DIR, "dge.chr19-diffs.local-genes.tsv"),
    sep = "\t",
    col.names = TRUE
)

fwrite(
    so_transcripts_cmpmt,
    file.path(DGE_DIR, "dge.chr19-diffs.local-transcripts.tsv"),
    sep = "\t",
    col.names = TRUE
)

saveRDS(
    so,
    file.path(DGE_DIR, "dge.chr19-diffs.sleuth-object.rds")
)
