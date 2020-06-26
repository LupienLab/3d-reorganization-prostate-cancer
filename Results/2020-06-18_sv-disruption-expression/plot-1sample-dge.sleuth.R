# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
suppressMessages(library("sleuth"))
source("../2020-02-19_chromoplexy/plotting-helper.R")

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
# load tests
sv_tests <- fread(
    "../2020-02-19_chromoplexy/Graphs/sv-disruption-tests.tsv",
    sep = "\t",
    header = TRUE
)
sv_tests$mut_samples <- split_comma_col(sv_tests$mut_samples)
sv_tests$nonmut_samples <- split_comma_col(sv_tests$nonmut_samples)

SAMPLES <- sv_tests[, sort(unique(unlist(mut_samples)))]

gencode <- fread("../../Data/External/GENCODE/gencode.v33.all-transcripts.bed", sep = "\t", header = FALSE, col.names = c("chr", "start", "end", "strand", "gene_id", "gene_name", "target_id", "transcript_name"))

# load all naive test results
tested_genes <- rbindlist(lapply(
    sv_tests$test_ID,
    function(tid) {
        dt <- fread(
            paste0("sleuth/test_", tid, ".genes.tested.tsv"),
            sep = "\t",
            header = TRUE
        )
        dt[, test_ID := tid]
        return(dt)
    }
))

tested_transcripts <- rbindlist(lapply(
    sv_tests$test_ID,
    function(tid) {
        dt <- fread(
            paste0("sleuth/test_", tid, ".transcripts.tested.tsv"),
            sep = "\t",
            header = TRUE
        )
        dt[, test_ID := tid]
        return(dt)
    }
))

transcript_table <- rbindlist(lapply(
    SAMPLES,
    function(s) {
        dt <- fread(
            paste0("../../Data/Processed/2020-06-17_PCa-RNA-seq/", s, "/abundance.tsv"),
            sep = "\t",
            header = TRUE
        )
        dt[, SampleID := s]
        return(dt)
    }
))


# ==============================================================================
# Analysis
# ==============================================================================
# combine transcripts for gene counts
gene_tpm <- rbindlist(apply(
    tested_genes,
    1,
    function(r) {
        # get gene_id
        gene <- r["target_id"]
        # get all transcripts associated with that gene
        transcript_ids <- gencode[gene_id == gene, target_id]
        # sum all TPM for all transcripts for that gene
        gene_tpm <- transcript_table[target_id %in% transcript_ids, .(tpm = sum(tpm)), by = "SampleID"]
        gene_tpm[, gene_id := gene]
        
        # get test_ID
        tid <- as.integer(r["test_ID"])
        # add Mutated column for each patient involved in a test
        gene_tpm[, test_ID := tid]
        gene_tpm[, Mutated := ""]
        mut_samples <- sv_tests[test_ID == tid, unlist(mut_samples)]
        nonmut_samples <- sv_tests[test_ID == tid, unlist(nonmut_samples)]
        gene_tpm[SampleID %in% mut_samples, Mutated := "Mutated"]
        gene_tpm[SampleID %in% nonmut_samples, Mutated := "Nonmutated"]

        # remove NA samples that are excluded from each test
        gene_tpm <- gene_tpm[Mutated != ""]
        return(gene_tpm)
    }
))

gene_tpm_mutated <- gene_tpm[, .(tpm = mean(tpm)), by = c("test_ID", "gene_id", "Mutated")]
gene_tpm_mutated_per_test <- dcast(gene_tpm_mutated, test_ID + gene_id ~ Mutated)
gene_tpm_mutated_per_test[, log2fc := log2((Mutated + 1) / (Nonmutated + 1))]

# merge with p-/q-values
merged_gene_tpm <- merge(
    x = tested_genes,
    y = gene_tpm_mutated_per_test,
    by.x = c("test_ID", "target_id"),
    by.y = c("test_ID", "gene_id")
)
merged_gene_tpm[, altered_group := (qval < 0.05 & abs(log2fc) > log(2))]
merged_gene_tpm[, test_ID := factor(test_ID)]

counted_gene_tpm <- rbindlist(lapply(
    sv_tests$test_ID,
    function(tid) {
        dt <- data.table(
            test_ID = tid,
            altered_group = c(FALSE, TRUE),
            N = c(
                merged_gene_tpm[test_ID == tid & altered_group == FALSE, .N],
                merged_gene_tpm[test_ID == tid & altered_group == TRUE, .N]
            ),
            Total = rep(merged_gene_tpm[test_ID == tid, .N], each = 2)
        )
        dt[, Frac := N / Total]
        return(dt)
    }
))

# ==============================================================================
# Plots
# ==============================================================================
gg_volcano <- (
    ggplot(data = tested_transcripts)
    + geom_point(aes(x = b * log2(exp(1)), y = -log10(qval), colour = (qval < 0.05 & abs(b) > log(2))))
    + labs(x = expression(log[2] * "(Fold change)"), y = expression(-log[10] * "(FDR)"))
    + scale_colour_manual(
        breaks = c(FALSE, TRUE),
        values = c("#BDBDBD", "#000000")
    )
    + theme_minimal()
)
savefig(gg_volcano, "Plots/volcano.transcripts")

gg_volcano <- (
    ggplot(data = merged_gene_tpm)
    + geom_point(aes(x = log2fc, y = -log10(qval), colour = (qval < 0.05 & abs(log2fc) > log2(2))))
    + geom_vline(aes(xintercept = -1))
    + geom_vline(aes(xintercept = 1))
    + labs(x = expression(log[2] * "(Fold change)"), y = expression(-log[10] * "(FDR)"))
    + scale_colour_manual(
        breaks = c(FALSE, TRUE),
        values = c("#BDBDBD", "#000000")
    )
    + coord_cartesian(
        xlim = c(-7, 7),
        ylim = c(0, 20)
    )
    + theme_minimal()
)
savefig(gg_volcano, "Plots/volcano.genes")

gg_test_exprs_dist <- (
    ggplot(data = counted_gene_tpm[complete.cases(counted_gene_tpm)])
    + geom_col(
        aes(
            x = factor(test_ID, ordered = TRUE, levels = counted_gene_tpm[altered_group == FALSE][order(Frac), unique(test_ID)]),
            y = Frac,
            fill = altered_group,
            group = factor(test_ID)
        ),
        position = "stack"
    )
    + labs(x = "Breakpoint", y = "Genes in TAD with breakpoint (%)")
    + scale_fill_discrete(
        breaks = c(FALSE, TRUE),
        labels = c("No difference", "Significant difference"),
        name = "Gene expression change"
    )
    + theme_minimal()
    + theme(
        axis.text.x = element_blank(),
        legend.position = "bottom"
    )
)
savefig(gg_test_exprs_dist, "Plots/barplot")
