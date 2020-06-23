# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("sleuth"))
suppressMessages(library("data.table"))
suppressMessages(library("argparse"))

if (!interactive()) {
    PARSER <- argparse::ArgumentParser(
        description = "Calculate gene expression fold changes using vanilla sleuth model"
    )
    PARSER$add_argument(
        "test_ID",
        type = "integer",
        help = "test_ID to plot comparison from"
    )
    ARGS <- PARSER$parse_args()
} else {
    ARGS <- list(
        "test_ID" = 0
    )
}

KALLISTO_DIR <- file.path("..", "..", "Data", "Processed", "2020-06-17_PCa-RNA-seq")


# ==============================================================================
# Functions
# ==============================================================================
split_comma_col <- function(v, f=identity) {
    # split into list
    splitv <- lapply(v, function(x) {strsplit(x, ",")[[1]]})
    # remove various non-informative characters (spaces, braces)
    splitv <- lapply(splitv, function(x) {gsub("[][ ]", "", x)})
    return(lapply(splitv, f))
}


# ==============================================================================
# Data
# ==============================================================================
# load sample metadata
metadata <- fread(
    "../../Data/External/LowC_Samples_Data_Available.tsv",
    sep = "\t",
    header = TRUE
)
metadata <- metadata[Include == "Yes"]
metadata[, SampleID :=  paste0("PCa", get("Sample ID"))]
SAMPLES <- metadata$SampleID

# add kallisto paths
metadata[, Path := file.path(KALLISTO_DIR, SampleID)]

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

# load TAD calls
tads <- fread(
    file.path("..", "2020-02-19_sv-disruption-TADs", "sv-disruption-tests.TADs.tsv"),
    sep = "\t",
    header = TRUE
)

# load test information
sv_tests <- fread(
    "../2020-02-19_chromoplexy/Graphs/sv-disruption-tests.tsv",
    sep = "\t",
    header = TRUE
)
sv_tests$mut_samples <- split_comma_col(sv_tests$mut_samples)
sv_tests$nonmut_samples <- split_comma_col(sv_tests$nonmut_samples)

# explicitly list mutated and non-mutated samples
mut_samples <- sv_tests[test_ID == ARGS$test_ID, unlist(mut_samples)]
nonmut_samples <- sv_tests[test_ID == ARGS$test_ID, unlist(nonmut_samples)]

# create experimental design for this test
design <- metadata[
    # only include samples involved in this test
    SampleID %in% c(mut_samples, nonmut_samples),
    # set default condition to be `Control`
    .(sample = SampleID, condition = "Control", path = Path)
]
# set condition of mut_samples to be `Mutated`
design[sample %in% mut_samples, condition := "Mutated"]


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
so <- sleuth_wt(so, "conditionMutated")

# save data
saveRDS(so, paste0("sleuth/test_", ARGS$test_ID, ".sleuth-object.rds"))

# extract results
so_genes <- as.data.table(sleuth_results(so, "conditionMutated", "wt", show_all = FALSE, pval_aggregate = TRUE))
# rows are duplicated with transcript IDs, but everything else is the same => unique only keeps relevant info
so_genes <- unique(so_genes[, .(target_id, gene_name, num_aggregated_transcripts, sum_mean_obs_counts, pval, qval)])

so_transcripts <- as.data.table(sleuth_results(so, "conditionMutated", "wt", show_all = FALSE, pval_aggregate = FALSE))
# NAs in start/end positions force the columns to be "character" instead of "integer"
# this happens because of the target mapping to GENCODE
# I've only included canonical chromosomes, not extra scaffolds, but the original index that kallisto quantifies against
# will contain these scaffolds, resulting in an NA target mapping from sleuth
so_transcripts <- so_transcripts[complete.cases(so_transcripts)]
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
    paste0("sleuth/test_", ARGS$test_ID, ".genes.all.tsv"),
    sep = "\t",
    col.names = TRUE
)
fwrite(
    so_transcripts,
    paste0("sleuth/test_", ARGS$test_ID, ".transcripts.all.tsv"),
    sep = "\t",
    col.names = TRUE
)

# focus on expression changes in TADs around region being tested
test_tad <- tads[test_ID == ARGS$test_ID]
test_genes <- so_genes[chr == test_tad$chr & start <= test_tad$end & end >= test_tad$start]
test_transcripts <- so_transcripts[chr == test_tad$chr & start <= test_tad$end & end >= test_tad$start]

fwrite(
    test_genes,
    paste0("sleuth/test_", ARGS$test_ID, ".genes.tested.tsv"),
    sep = "\t",
    col.names = TRUE
)
fwrite(
    test_transcripts,
    paste0("sleuth/test_", ARGS$test_ID, ".transcripts.tested.tsv"),
    sep = "\t",
    col.names = TRUE
)
