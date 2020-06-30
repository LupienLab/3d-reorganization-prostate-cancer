# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("stringr"))
suppressMessages(library("ggplot2"))
suppressMessages(library("argparse"))
source("../2020-02-19_chromoplexy/plotting-helper.R")

if (!interactive()) {
    PARSER <- argparse::ArgumentParser(
        description = "Plot the TPM comparison for a given gene in a given test"
    )
    PARSER$add_argument(
        "test_ID",
        type = "integer",
        help = "test_ID to plot comparison from"
    )
    PARSER$add_argument(
        "loc",
        type = "character",
        help = "Genes overlapping this position will be plotted",
    )
    ARGS <- PARSER$parse_args()
} else {
    ARGS <- list(
        "test_ID" = 44,
        "loc" = "chr14:35340000-35930000"
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
cat("Loading data\n")
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

# full table of estimated counts as well as T2E-comparison results
so_transcripts <- fread("results.transcripts.tsv", sep = "\t")
so_genes <- fread("results.genes.tsv", sep = "\t")
full_table <- fread("all-samples.abundance.tsv", sep = "\t")

# GENCODE transcript annotation
gencode <- fread(
    file.path("..", "..", "Data", "External", "GENCODE", "gencode.v33.all-transcripts.bed"),
    sep = "\t",
    header = FALSE,
    col.names = c("chr", "start", "end", "strand", "gene_id", "gene_name", "transcript_id", "transcript_name")
)

# ==============================================================================
# Analysis
# ==============================================================================
cat("Getting gene information\n")
pos_regex <- str_match(ARGS$loc, "^(chr[0-9]{0,3}[XYM]?):(\\d+)-(\\d+)$")
pos <- list(
    chr = pos_regex[, 2],
    start = as.integer(pos_regex[, 3]),
    end = as.integer(pos_regex[, 4])
)
# get transcript IDs overlapping this posiiton
overlapping_gencode <- gencode[(chr == pos$chr) & (end >= pos$start) & (start <= pos$end)]

# get query transcript counts in TPM
query_transcript_tpm <- rbindlist(mapply(
    function(target, gene_id, gene_name) {
        dt <- full_table[target_id == target, .SD]
        dt[, gene_id := gene_id]
        dt[, gene_name := gene_name]
        return(dt)
    },
    overlapping_gencode[, transcript_id],
    overlapping_gencode[, gene_id],
    overlapping_gencode[, gene_name],
    SIMPLIFY = FALSE
))
colnames(query_transcript_tpm)[colnames(query_transcript_tpm) == "sample"] <- "SampleID"
query_transcript_tpm[, log10p1TPM := log10(tpm + 1)]

# get query gene counts in TPM
query_gene_tpm <- query_transcript_tpm[, .(tpm = sum(tpm)), by = c("SampleID", "gene_id", "gene_name")]

this_test_stratification <- data.table(
    SampleID = c(mut_samples, nonmut_samples),
    Mutated = rep(c("Yes", "No"), c(length(mut_samples), length(nonmut_samples)))
)

# merge test group information for the desired test
query_gene_tpm <- merge(
    x = query_gene_tpm,
    y = this_test_stratification
)
query_gene_tpm[, log10p1TPM := log10(tpm + 1)]

# ==============================================================================
# Plots
# ==============================================================================
cat("Plotting\n")
# merge sample colour information
query_gene_tpm <- merge(
    x = query_gene_tpm,
    y = metadata[, .SD, .SDcols = c("SampleID", "Colour")]
)
query_transcript_tpm <- merge(
    x = query_gene_tpm,
    y = metadata[, .SD, .SDcols = c("SampleID", "Colour")]
)

paths <- rbindlist(lapply(
    query_gene_tpm[, unique(gene_name)],
    function(gene) {
        dt <- data.table(
            gene_name = gene,
            Mutated = rep(c("No", "Yes"), each = 2),
            y = c(
                (query_gene_tpm[Mutated == "No" & gene_name == gene, max(log10p1TPM)] + query_gene_tpm[gene_name == gene, max(log10p1TPM) * 1.05]) / 2,
                query_gene_tpm[gene_name == gene, max(log10p1TPM) * 1.05],
                query_gene_tpm[gene_name == gene, max(log10p1TPM) * 1.05],
                (query_gene_tpm[Mutated == "Yes" & gene_name == gene, max(log10p1TPM)] + query_gene_tpm[gene_name == gene, max(log10p1TPM) * 1.05]) / 2
            )
        )
    }
))

gg_gene_tpm <- (
    ggplot(data = query_gene_tpm)
    + geom_boxplot(
        aes(x = gene_name, y = log10p1TPM, fill = Mutated, group = paste(Mutated, gene_name, sep = "_")),
        alpha = 0.3,
        colour = "black",
        outlier.shape = NA
    )
    + geom_point(
        aes(x = gene_name, y = log10p1TPM, colour = Colour, group = paste(Mutated, gene_name, sep = "_")),
        position = position_dodge(width = 0.8)
    )
    # + geom_path(
    #     data = paths,
    #     mapping = aes(x = gene_name, y = y, group = paste(Mutated, gene_name, sep = "_")),
    #     colour = "black"
    # )
    + scale_x_discrete(
        name = "Gene"
    )
    + scale_fill_manual(
        breaks = c("No", "yes"),
        values = c("#ffffff", "#424242"),
        name = NULL
    )
    + scale_colour_manual(
        breaks = query_gene_tpm[, unique(Colour)],
        values = query_gene_tpm[, unique(Colour)],
        name = NULL
    )
    + labs(y = expression(log[10] * " (TPM + 1)"))
    + guides(colour = FALSE, fill = FALSE)
    + coord_flip()
    + theme_minimal()
    + theme(
    )
)
savefig(gg_gene_tpm, paste0("Plots/", paste(pos, collapse = "_"), ".expression.gene-level"), height = 10, width = 8)
