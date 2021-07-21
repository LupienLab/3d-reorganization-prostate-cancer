# ==============================================================================
# Meta
# ==============================================================================
# Plot gene expression
# --------------------------------------
# Description: Plot the gene expression changes for the T2E translocation insertion site
# Author: James Hawley


# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("logging"))

loginfo("Loading packages")
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
suppressMessages(library("sleuth"))
source(file.path("..", "src", "savefig.R"))
source("helper-functions.R")

RES_DIR <- file.path(
    "..", "..", "results", "2020-06-18_sv-disruption-expression"
)
PLOT_DIR <- file.path(RES_DIR, "Plots")


# ==============================================================================
# Data
# ==============================================================================
loginfo("Loading data")

# mutated/non-mutated samples
mut_sample <- "PCa53687"
nonmut_samples <- setdiff(
    c(
        "PCa13266", "PCa13848", "PCa14121", "PCa19121", "PCa3023", "PCa33173",
        "PCa40507", "PCa51852", "PCa53687", "PCa56413", "PCa57294", "PCa58215"
    ),
    mut_sample
)

# load GENCODE annotations
gencode_genes <- fread(
    file.path("..", "..", "data", "External", "GENCODE", "gencode.v33.all-genes.bed"),
    sep = "\t",
    header = FALSE,
    col.names = c("chr", "start", "end", "strand", "ens_gene", "gene_name")
)
gencode_tx <- fread(
    file.path("..", "..", "data", "External", "GENCODE", "gencode.v33.all-transcripts.bed"),
    sep = "\t",
    header = FALSE,
    col.names = c("chr", "start", "end", "strand", "ens_gene", "gene_name", "ens_transcript", "transcript_name")
)

# ==============================================================================
# Analysis
# ==============================================================================
loginfo("Performing calculations")

# load Kallisto/sleuth counts
so_counts <- get_important_sleuth_info(
    path = file.path(RES_DIR, "sleuth", "test_201.sleuth-object.rds"),
    alt_gene = "BRAF",
    base = 2
)

tx_counts <- so_counts$tx_counts

# set order of transcripts by mean expression of non-mutant samples
tx_order <- tx_counts[
    sample != mut_sample,
    .(mean_est_counts = mean(est_counts)),
    by = target_id
][,
    target_id,
    keyby = "mean_est_counts"
]$target_id
tx_counts[, target_id := factor(target_id, ordered = TRUE, levels = rev(tx_order))]


# ==============================================================================
# Plots
# ==============================================================================
loginfo("Plotting data")

gg <- (
    ggplot()
    + geom_boxplot(
        data = tx_counts[sample != mut_sample],
        mapping = aes(
            x = target_id,
            y = est_counts
        ),
        outlier.shape = NA,
        alpha = 0.7
    )
    + geom_point(
        data = tx_counts[sample != mut_sample],
        mapping = aes(
            x = target_id,
            y = est_counts
        ),
        fill = "#BDBDBD",
        colour = "black",
        position = position_jitter(height = 0, width = 0.3),
        shape = 21,
        size = 4,
        alpha = 0.5
    )
    + geom_point(
        data = tx_counts[sample == mut_sample],
        mapping = aes(
            x = target_id,
            y = est_counts
        ),
        fill = "#ff6347",
        colour = "black",
        shape = 21,
        size = 4
    )
    + labs(
        x = "BRAF Transcript",
        y = expression("Estimated mRNA count (" * log[2] * ")")
    )
    + theme_minimal()
    + jrh_theme()
)
savefig(gg, file.path(PLOT_DIR, "tx-exprs.BRAF"), height = 10, width = 6)

# ==============================================================================
# Save data
# ==============================================================================
# save the BRAF transcript data
fwrite(
    tx_counts,
    file.path(RES_DIR, "quant.BRAF.tsv"),
    sep = "\t",
    col.names = TRUE
)

# save differential analysis data
fwrite(
    so_counts$tx_de,
    file.path(RES_DIR, "test.BRAF.tsv"),
    sep = "\t",
    col.names = TRUE
)
