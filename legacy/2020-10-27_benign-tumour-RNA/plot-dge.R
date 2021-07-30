# ==============================================================================
# Meta
# ==============================================================================
# plot-dge
# --------------------------------------
# Description: tumour-vs-benign
# Author: James Hawley

suppressMessages(library("sleuth"))
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
source(file.path("..", "src", "savefig.R"))


# ==============================================================================
# Data
# ==============================================================================
# load sample metadata
metadata <- fread(
    "config.tsv",
    sep = "\t",
    header = TRUE
)
metadata <- metadata[Include == "Yes"]
SAMPLES <- metadata$SampleID

so_transcripts <- fread("results.transcripts.tsv", sep = "\t")
so_genes <- fread("results.genes.tsv", sep = "\t")

transcripts <- fread(
    file.path("..", "..", "Data", "External", "GENCODE", "gencode.v33.all-transcripts.bed"),
    sep = "\t",
    header = FALSE,
    col.names = c("chr", "start", "end", "strand", "ens_gene", "gene_name", "target_id", "transcript_name")
)

# find genes with at least one transcript with a > 2 x fold change
# log(2) since b coefficients are on natural log estimates of counts
large_fc_tx <- so_transcripts[qval < 0.01 & abs(b) > log(2), target_id]
large_fc_genes <- transcripts[target_id %in% large_fc_tx, ens_gene]

large_fc_sig_genes <- so_genes[target_id %in% large_fc_genes & qval < 0.01]
fwrite(
    large_fc_sig_genes[, .SD, .SDcols = c(
        "chr", "start", "end", "gene_name", "qval", "strand",
        "target_id", "num_aggregated_transcripts", "sum_mean_obs_counts", "pval"
    )],
    "results.genes.sig.bed",
    sep = "\t",
    col.names = TRUE
)

# ==============================================================================
# Plots
# ==============================================================================
gg_volcano <- (
    ggplot(data = so_transcripts)
    +
        geom_point(
            aes(
                x = b / log(2),
                y = -log10(pval),
                fill = paste(
                    qval < 0.05,
                    abs(b / log(2)) > 1,
                    sign(b)
                )
            ),
            shape = 21,
            alpha = 0.01
        )
        + scale_x_continuous(
            name = bquote(log[2] * "(Expression Fold Change) (Tumour / Benign)")
        )
        + scale_y_continuous(
            name = bquote(-log[10] * "(p-value)")
        )
        + scale_fill_manual(
            name = NULL,
            breaks = c(
                "FALSE FALSE -1",
                "FALSE FALSE 1",
                "FALSE TRUE -1",
                "FALSE TRUE 1",
                "TRUE FALSE -1",
                "TRUE FALSE 1",
                "TRUE TRUE -1",
                "TRUE TRUE 1"
            ),
            values = rep(
                c("#BDBDBD", "#6495ED", "#FFA500"),
                c(6, 1, 1)
            )
        )
        + guides(colour = FALSE, fill = FALSE)
        # + facet_wrap(~ target_id, drop = TRUE, scale = "free")
        +
        theme_minimal()
    # + theme(
    #     axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 1)
    # )
)
savefig(gg_volcano, "Plots/volcano.transcripts")

gg_pval_hist_genes <- (
    ggplot(
        data = so_genes,
        mapping = aes(x = pval)
    )
    + geom_histogram()
    + theme_minimal()
)
savefig(
    gg_pval_hist_genes,
    file.path(PLOT_DIR, "pvals.genes")
)
gg_pval_hist_tx <- (
    ggplot(
        data = so_transcripts,
        mapping = aes(x = pval)
    )
    + geom_histogram()
    + theme_minimal()
)
savefig(
    gg_pval_hist_tx,
    file.path(PLOT_DIR, "pvals.transcripts")
)

gg_erg_sum <- (
    ggplot(data = erg_tpm)
    +
        geom_boxplot(
            aes(x = condition, y = sum(TPM), fill = condition),
            alpha = 0.3,
            outlier.shape = NA
        )
        +
        geom_point(
            aes(x = condition, y = sum(TPM), colour = condition),
            position = position_jitter(height = 0, width = 0.2)
        )
        +
        scale_x_discrete(
            breaks = c("No", "Yes"),
            labels = c("T2E-", "T2E+"),
            name = NULL
        )
        +
        guides(colour = FALSE, fill = FALSE)
        +
        theme_minimal()
        +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 1)
        )
)
savefig(gg_erg_sum, "Plots/ERG-expression.gene-level")

gg_znf_sum <- (
    ggplot(data = znf_tpm[, .(TPM = sum(TPM)), by = c("condition", "SampleID")])
    +
        geom_boxplot(
            aes(x = condition, y = TPM, fill = condition),
            alpha = 0.3,
            outlier.shape = NA
        )
        +
        geom_point(
            aes(x = condition, y = TPM, colour = condition),
            position = position_jitter(height = 0, width = 0.2)
        )
        +
        geom_path(
            data = data.table(
                x = rep(c("No", "Yes"), each = 2),
                y = c(10, 25, 25, 24.5),
                group = 1
            ),
            mapping = aes(x = x, y = y, group = group)
        )
        +
        geom_text(
            data = data.table(
                x = 1.5,
                y = 25,
                label = paste0("p = ", so_genes[gene_name == "ZNF516", signif(pval, digits = 3)])
            ),
            mapping = aes(x = x, y = y, label = label),
            vjust = -0.8
        )
        +
        scale_x_discrete(
            breaks = c("No", "Yes"),
            labels = c("T2E-", "T2E+"),
            name = NULL
        )
        +
        scale_colour_manual(
            breaks = c("No", "Yes"),
            labels = c("T2E-", "T2E+"),
            values = c("#418B3D", "#3215C1"),
            name = NULL
        )
        +
        guides(colour = FALSE, fill = FALSE)
        +
        theme_minimal()
        +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 1)
        )
)
savefig(gg_znf_sum, "Plots/ZNF516-expression.gene-level", width = 8)