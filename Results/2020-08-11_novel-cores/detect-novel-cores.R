# ==============================================================================
# Meta
# ==============================================================================
# Detect novel CORES
# --------------------------------------
# Description: Find instances of novel COREs around SV breakpoints
# Author: James Hawley


# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("regioneR"))

GRAPH_DIR <- file.path("..", "2020-02-19_chromoplexy", "Graphs")
CORE_DIR <- file.path("..", "..", "Data", "Processed", "2019-05-03_PCa-H3K27ac-peaks", "COREs")
EXPRS_DIR <- file.path("..", "2020-06-18_sv-disruption-expression", "sleuth")
ACETYL_DIR <- file.path("..", "2020-06-12_sv-disruption-acetylation", "Acetylation", "Tests")

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
metadata <- fread("config.tsv", sep = "\t", header = TRUE)
SAMPLES <- metadata[, SampleID]

# load differential analysis results
res <- fread(
    file.path("..", "2020-06-12_sv-disruption-acetylation", "exprs-acetyl.tsv"),
    sep = "\t",
    header = TRUE
)

# load SV test information
tests <- fread(file.path(GRAPH_DIR, "sv-disruption-tests.tsv"), sep = "\t", header = TRUE)
tests$mut_samples <- split_comma_col(tests$mut_samples)
tests$nonmut_samples <- split_comma_col(tests$nonmut_samples)

# load gene locations
gencode <- fread(
    file.path("..", "..", "Data", "External", "GENCODE", "gencode.v33.all-genes.bed"),
    sep = "\t",
    header = FALSE,
    col.names = c("chr", "start", "end", "strand", "ens_gene", "gene_name")
)

# load sample COREs
cores <- lapply(
    SAMPLES,
    function(s) {
        dt <- fread(file.path(CORE_DIR, paste0(s, ".cores.bed")), sep = "\t", header = FALSE, col.names = c("chr", "start", "end"))
        gr <- toGRanges(dt, genome = "hg38")
        return(gr)
    }
)
names(cores) <- SAMPLES


# ==============================================================================
# Analysis
# ==============================================================================
res_pos <- merge(
    x = gencode,
    y = res,
    by = c("ens_gene", "gene_name")
)
setcolorder(res_pos, c("chr", "start", "end"))

res_gr <- toGRanges(res_pos, genome = "hg38")
strand(res_gr) <- res_pos$strand
mcols(res_gr)$strand <- NULL

counters <- list(
    "COREs_near_breakpoint" = 0,
    "COREs_near_breakpoint_test_IDs" = c(),
    "COREs_near_breakpoint_genes" = c(),
    "COREs_near_breakpoint_log2exprs_fc" = c(),
    "novel_COREs" = 0,
    "novel_CORE_test_IDs" = c(),
    "novel_CORE_genes" = c(),
    "novel_CORE_log2exprs_fc" = c()
)

# iterate over differentially expressed genes
for (i in 1:length(res_gr)) {
    row <- res_gr[i]
    if (row$exprs_qval >= 0.05) next
    ti <- row$test_ID
    # get (non-)mutated samples from this test
    mut_samples <- tests[test_ID == ti, unlist(mut_samples)]
    nonmut_samples <- tests[test_ID == ti, unlist(nonmut_samples)]
    # get COREs from mutated and non-mutated samples
    mut_cores <- reduce(unlist(GRangesList(lapply(mut_samples, function(s) cores[[s]]))))
    nonmut_cores <- reduce(unlist(GRangesList(lapply(nonmut_samples, function(s) cores[[s]]))))

    # check if gene in question overlaps a CORE from this sample
    deg_core_overlap <- findOverlaps(row, mut_cores)
    if (length(deg_core_overlap) == 0) next
    
    # record that there is a CORE near the breakpoint
    counters[["COREs_near_breakpoint"]] <- counters[["COREs_near_breakpoint"]] + 1
    counters[["COREs_near_breakpoint_test_IDs"]] <- c(counters[["COREs_near_breakpoint_test_IDs"]], ti)
    counters[["COREs_near_breakpoint_genes"]] <- c(counters[["COREs_near_breakpoint_genes"]], row$gene_name)
    counters[["COREs_near_breakpoint_log2exprs_fc"]] <- c(counters[["COREs_near_breakpoint_log2exprs_fc"]], row$log2exprs_fc)
    
    # check that the CORE is novel (not found in any of the non-mutated samples)
    deg_nonmut_core_overlap <- findOverlaps(row, nonmut_cores)

    if (length(deg_nonmut_core_overlap) == 0) {
        counters[["novel_COREs"]] <- counters[["novel_COREs"]] + 1
        counters[["novel_CORE_test_IDs"]] <- c(counters[["novel_CORE_test_IDs"]], ti)
        counters[["novel_CORE_genes"]] <- c(counters[["novel_CORE_genes"]], row$gene_name)
        counters[["novel_CORE_log2exprs_fc"]] <- c(counters[["novel_CORE_log2exprs_fc"]], row$log2exprs_fc)
    }
    # cat("Test ID: ", ti, "\n")
    # cat("Gene:\n")
    # print(row)
    # cat("\nMutated CORE:\n")
    # print(mut_cores[subjectHits(deg_core_overlap)])
    # cat("\nNon-mutated CORE:\n")
    # print(nonmut_cores[subjectHits(deg_nonmut_core_overlap)])
}
counters
