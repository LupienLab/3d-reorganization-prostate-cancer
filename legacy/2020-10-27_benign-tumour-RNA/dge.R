# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("sleuth"))
suppressMessages(library("logging"))

KALLISTO_DIR <- file.path("..", "..", "Data", "Processed", "2020-06-17_PCa-RNA-seq")

# ==============================================================================
# Data
# ==============================================================================
loginfo("Loading data")
# load sample metadata
metadata <- fread(
	"config.tsv",
	sep = "\t",
	header = TRUE
)
metadata <- metadata[Include == "Yes"]
SAMPLES <- metadata$SampleID

# add kallisto paths
metadata[, Path := file.path(KALLISTO_DIR, SampleID)]

design <- metadata[, .(sample = SampleID, condition = Type, path = Path)]

# load GENCODE annotations
genes <- fread(
	file.path("..", "..", "Data", "External", "GENCODE", "gencode.v33.all-genes.bed"),
	sep = "\t",
	header = FALSE,
	col.names = c("chr", "start", "end", "strand", "ens_gene", "gene_name")
)
transcripts <- fread(
	file.path("..", "..", "Data", "External", "GENCODE", "gencode.v33.all-transcripts.bed"),
	sep = "\t",
	header = FALSE,
	col.names = c("chr", "start", "end", "strand", "ens_gene", "gene_name", "target_id", "transcript_name")
)


# ==============================================================================
# Analysis
# ==============================================================================
loginfo("Preparing sleuth object")
# create sleuth object, aggregating transcripts to genes
so <- sleuth_prep(
	design,
	extra_bootstrap_summary = TRUE,
	target_mapping = transcripts,
	aggregation_column = "ens_gene",
	num_cores = 2
)

loginfo("Fitting the model")
# fit full model
so <- sleuth_fit(so, ~condition, "full")

# fit reduced (null) model
so <- sleuth_fit(so, ~1, "reduced")

# perform differential analysis
so <- sleuth_wt(so, "conditionMalignant")

# save data
saveRDS(so, "sleuth-object.rds")

# extract results
so_genes <- as.data.table(sleuth_results(so, "conditionMalignant", "wt", show_all = FALSE, pval_aggregate = TRUE))
# rows are duplicated with transcript IDs, but everything else is the same => unique only keeps relevant info
so_genes <- unique(so_genes[, .(target_id, gene_name, num_aggregated_transcripts, sum_mean_obs_counts, pval, qval)])

so_transcripts <- as.data.table(sleuth_results(so, "conditionMalignant", "wt", show_all = FALSE, pval_aggregate = FALSE))
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
fwrite(so_genes, "results.genes.tsv", sep = "\t", col.names = TRUE)
fwrite(so_transcripts, "results.transcripts.tsv", sep = "\t", col.names = TRUE)

# write full sleuth count table
full_table <- as.data.table(kallisto_table(so))
fwrite(
	full_table[, .SD, .SDcols = c("target_id", "sample", "est_counts", "tpm", "eff_len", "len")],
	"all-samples.abundance.tsv",
	sep = "\t",
	col.names = TRUE
)

