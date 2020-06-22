# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("sleuth"))
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
suppressMessages(library("argparse"))
source("../2020-02-19_chromoplexy/plotting-helper.R")

if (!interactive()) {
    PARSER <- argparse::ArgumentParser(
        description = "Plot the TPM comparison for a given gene in a given test"
    )
    PARSER$add_argument(
        "gene",
        type = "character",
        help = "Gene name to be shown"
    )
    PARSER$add_argument(
        "test_ID",
        type = "integer",
        help = "test_ID to plot comparison from"
    )
    ARGS <- PARSER$parse_args()
} else {
    ARGS <- list(
        "gene" = "ERG",
        "test_ID" = 46
    )
}

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

# load test information
tests <- fread(
    "../2020-02-19_chromoplexy/Graphs/sv-disruption-tests.tsv",
    sep = "\t",
    header = TRUE
)
tests$mut_samples <- split_comma_col(tests$mut_samples)
tests$nonmut_samples <- split_comma_col(tests$nonmut_samples)

# explicitly list mutated and non-mutated samples
mut_samples <- tests[test_ID == ARGS$test_ID, unlist(mut_samples)]
nonmut_samples <- tests[test_ID == ARGS$test_ID, unlist(nonmut_samples)]

so <- readRDS("sleuth-object.rds")
so_transcripts <- fread("results.transcripts.tsv", sep = "\t")
so_genes <- fread("results.genes.tsv", sep = "\t")

# full table of estimated counts
full_table <- as.data.table(kallisto_table(so))

# ==============================================================================
# Analysis
# ==============================================================================
query_transcript_counts <- rbindlist(lapply(
    so_transcripts[gene_name == ARGS$gene, target_id],
    function(target) {
        dt <- as.data.table(get_bootstrap_summary(so, target, "est_counts"))
        dt[, target_id := target]
        return(dt)
    }
))

# get query transcript counts in TPM
query_transcript_tpm <- merge(
    x = query_transcript_counts,
    y = full_table[, .(target_id, sample, tpm)],
    by = c("target_id", "sample")
)
colnames(query_transcript_tpm)[colnames(query_transcript_tpm) == "sample"] <- "SampleID"
colnames(query_transcript_tpm)[colnames(query_transcript_tpm) == "tpm"] <- "TPM"
query_transcript_tpm[, logTPM := log10(TPM + 1)]

# get query gene counts in TPM
query_gene_tpm <- query_transcript_tpm[, .(TPM = sum(TPM)), by = "SampleID"]


this_test_stratification <- data.table(
    SampleID = c(mut_samples, nonmut_samples),
    Mutated = rep(c("Yes", "No"), c(length(mut_samples), length(nonmut_samples)))
)

# merge test group information for the desired test
query_gene_tpm <- merge(
    x = query_gene_tpm,
    y = this_test_stratification
)
query_gene_tpm[, logTPM := log10(TPM + 1)]

# ==============================================================================
# Plots
# ==============================================================================
# merge sample colour information
query_gene_tpm <- merge(
    x = query_gene_tpm,
    y = metadata[, .SD, .SDcols = c("SampleID", "Colour")]
)

gg_gene_tpm <- (
    ggplot(data = query_gene_tpm)
    + geom_boxplot(
        aes(x = Mutated, y = logTPM, fill = Mutated),
        alpha = 0.3,
        colour = "black",
        outlier.shape = NA
    )
    + geom_point(aes(x = Mutated, y = logTPM, colour = Colour), position = position_jitter(height = 0, width = 0.2))
    + geom_path(
        data = data.table(
            x = rep(c("No", "Yes"), each = 2),
            y = c(
                (query_gene_tpm[Mutated == "No", max(logTPM)] + query_gene_tpm[, max(logTPM) * 1.05]) / 2,
                query_gene_tpm[, max(logTPM) * 1.05],
                query_gene_tpm[, max(logTPM) * 1.05],
                (query_gene_tpm[Mutated == "Yes", max(logTPM)] + query_gene_tpm[, max(logTPM) * 1.05]) / 2
            ),
            group = 1
        ),
        mapping = aes(x = x, y = y, group = group),
        colour = "black"
    )
    # + geom_text(
    #     data = data.table(
    #         x = 1.5,
    #         y = query_gene_tpm[, max(TPM) + 1],
    #         label = paste0("p = ", so_genes[gene_name == ARGS$gene, signif(pval, digits = 3)])
    #     ),
    #     mapping = aes(x = x, y = y, label = label),
    #     vjust = -0.8
    # )
    + scale_x_discrete(
        breaks = c("No", "Yes"),
        name = "Mutated"
    )
    + scale_fill_manual(
        breaks = c("No", "yes"),
        values = c("#ffffff", "#424242"),
        name = NULL
    )
    + scale_colour_manual(
        breaks = query_gene_tpm$Colour,
        values = query_gene_tpm$Colour,
        name = NULL
    )
    + labs(
        y = expression(log[10] * " (TPM + 1)"),
        title = ARGS$gene
    )
    + guides(colour = FALSE, fill = FALSE)
    + theme_minimal()
    + theme(
        axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 1)
    )
)
savefig(gg_gene_tpm, paste0("Plots/", ARGS$gene, "-expression.gene-level"), width = 8)
