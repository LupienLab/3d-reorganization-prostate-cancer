# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("sleuth"))
suppressMessages(library("data.table"))
source(
    file.path("..", "2020-02-19_chromoplexy", "plotting-helper.R")
)

KALLISTO_DIR <- file.path("..", "..", "data", "Processed", "2020-06-17_PCa-RNA-seq")
RES_DIR <- file.path(
    "..", "..", "results", "2020-06-18_sv-disruption-expression"
)
PLOT_DIR <- file.path(RES_DIR, "Plots")

# ==============================================================================
# Data
# ==============================================================================
# load sample metadata
metadata <- fread(
    file.path(
        "..", "..", "data", "External", "LowC_Samples_Data_Available.tsv"
    ),
    sep = "\t",
    header = TRUE
)
metadata <- metadata[Include == "Yes"]
metadata[, SampleID :=  paste0("PCa", get("Sample ID"))]
SAMPLES <- metadata$SampleID

# add kallisto paths
metadata[, Path := file.path(KALLISTO_DIR, SampleID)]

design <- metadata[, .(sample = SampleID, condition = get("T2E Status"), path = Path)]

# load GENCODE annotations
genes <- fread(
    file.path("..", "..", "data", "External", "GENCODE", "gencode.v33.all-genes.bed"),
    sep = "\t",
    header = FALSE,
    col.names = c("chr", "start", "end", "strand", "ens_gene", "gene_name")
)
transcripts <- fread(
    file.path("..", "..", "data", "External", "GENCODE", "gencode.v33.all-transcripts.bed"),
    sep = "\t",
    header = FALSE,
    col.names = c("chr", "start", "end", "strand", "ens_gene", "gene_name", "target_id", "transcript_name")
)

tads <- fread(
    file.path("..", "..", "results", "2020-02-19_sv-disruption-TADs", "sv-disruption-tests.TADs.tsv"),
    sep = "\t",
    header = TRUE
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
    num_cores = 2
)

# fit full model
so <- sleuth_fit(so, ~condition, "full")

# fit reduced (null) model
so <- sleuth_fit(so, ~1, "reduced")

# perform differential analysis
so <- sleuth_lrt(so, "reduced", "full")
so <- sleuth_wt(so, "conditionYes")

# save data
saveRDS(so, file.path(RES_DIR, "sleuth-object.rds"))

# extract results
so_genes <- as.data.table(sleuth_results(so, "conditionYes", "wt", show_all = FALSE, pval_aggregate = TRUE))
# rows are duplicated with transcript IDs, but everything else is the same => unique only keeps relevant info
so_genes <- unique(so_genes[, .(target_id, gene_name, num_aggregated_transcripts, sum_mean_obs_counts, pval, qval)])

so_transcripts <- as.data.table(sleuth_results(so, "conditionYes", "wt", show_all = FALSE, pval_aggregate = FALSE))
# NAs in start/end positions force the columns to be "character" instead of "integer"
# this happens because of the target mapping to GENCODE
# I've only included canonical chromosomes, not extra scaffolds, but the original index that kallisto quantifies against
# will contain these scaffolds, resulting in an NA target mapping from sleuth
so_transcripts[, `:=`(start = as.integer(start), end = as.integer(end))]


# merge gene-level annotations
so_genes <- merge(
    x = so_genes,
    y = genes,
    by.x = c("target_id", "gene_name"),
    by.y = c("ens_gene", "gene_name"),
    all.x = TRUE
)[order(qval)]

# save annotated table
fwrite(
    so_genes,
    file.path(RES_DIR, "results.genes.tsv"),
    sep = "\t",
    col.names = TRUE
)
fwrite(
    so_transcripts,
    file.path(RES_DIR, "results.transcripts.tsv"),
    sep = "\t",
    col.names = TRUE
)

# sleuth_live(so, options = list(port = 42427, launch.browser = FALSE))

# compare expression changes in TADs around either breakpoint for the T2E fusion
# endpoints for T2E are in test_ID 45 (TMPRSS2) and 46 (ERG)
erg_tad <- tads[test_ID == 46]
tmprss2_tad <- tads[test_ID == 45]
t2e_genes <- rbindlist(list(
    # get transcripts from within ERG's TAD
    so_genes[chr == erg_tad$chr & start <= erg_tad$end & end >= erg_tad$start],
    so_genes[chr == tmprss2_tad$chr & start <= tmprss2_tad$end & end >= tmprss2_tad$start]
))
t2e_transcripts <- rbindlist(list(
    # get transcripts from within ERG's TAD
    so_transcripts[chr == erg_tad$chr & start <= erg_tad$end & end >= erg_tad$start],
    so_transcripts[chr == tmprss2_tad$chr & start <= tmprss2_tad$end & end >= tmprss2_tad$start]
))

fwrite(
    t2e_genes,
    file.path(RES_DIR, "t2e-comparison.genes.tsv"),
    sep = "\t",
    col.names = TRUE
)
fwrite(
    t2e_transcripts,
    file.path(RES_DIR, "t2e-comparison.transcripts.tsv"),
    sep = "\t",
    col.names = TRUE
)

# write full sleuth count table
full_table <- as.data.table(kallisto_table(so))
fwrite(
    full_table[, .SD, .SDcols = c("target_id", "sample", "est_counts", "tpm", "eff_len", "len")],
    file.path(RES_DIR, "all-samples.abundance.tsv"),
    sep = "\t",
    col.names = TRUE
)
