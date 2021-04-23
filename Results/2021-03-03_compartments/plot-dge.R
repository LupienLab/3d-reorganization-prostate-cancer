# ==============================================================================
# Meta
# ==============================================================================
# plot-dge
# ------------------------------------------------
# Author: James Hawley
# Description: Plot differential gene expression results in regions with differential compartments


# ==============================================================================
# Environment
# ==============================================================================
suppressWarnings(library("logging"))
loginfo("Loading packages")

suppressWarnings(library("data.table"))
suppressWarnings(library("ggplot2"))
suppressWarnings(library("ggrepel"))

CMPMT_DIR <- "Compartments"
DGE_DIR <- "DGE"
PLOT_DIR <- "Plots"


# ==============================================================================
# Functions
# ==============================================================================
# Helper function for plotting
labeller <- function(s) {
    ifelse(
        s == "benign",
        "Benign-Like",
        ifelse(
            s == "novel",
            "Altered",
            s
        )
    )
}

# ==============================================================================
# Data
# ==============================================================================
loginfo("Loading data")

# load sample metadata
meta <- fread("config.tsv", sep = "\t", header = TRUE)
meta <- meta[Include == "Yes"]
SAMPLES <- meta[Type == "Malignant", Sample_ID]

# load chr19 differential expression
chr19_genes <- fread(
    file.path(DGE_DIR, "dge.chr19-diffs.local-genes.tsv"),
    sep = "\t",
    header = TRUE
)
chr19_tx <- fread(
    file.path(DGE_DIR, "dge.chr19-diffs.local-transcripts.tsv"),
    sep = "\t",
    header = TRUE
)
so_chr19 <- readRDS(file.path(DGE_DIR, "dge.chr19-diffs.sleuth-object.rds"))

dge_gene_ids <- chr19_genes[qval < 0.05, target_id]

# load chrY differential expression
chrY_genes <- fread(
    file.path(DGE_DIR, "dge.chrY-diffs.local-genes.tsv"),
    sep = "\t",
    header = TRUE
)
chrY_tx <- fread(
    file.path(DGE_DIR, "dge.chrY-diffs.local-transcripts.tsv"),
    sep = "\t",
    header = TRUE
)

# ==============================================================================
# Plots
# ==============================================================================
loginfo("Plotting figures")

for (gene_id in dge_gene_ids) {
    # get the transcript IDs for this gene, sorted by qval for easy plotting later
    tx_ids <- chr19_tx[ens_gene == gene_id, target_id, keyby = "qval"]
    gene_name <- chr19_genes[target_id == gene_id, gene_name]
    df <- rbindlist(lapply(
        tx_ids[, target_id],
        function(tid) {
            df1 <- as.data.table(get_bootstrap_summary(so_chr19, tid))
            df1[, transcript_id := tid]
            df1[, transcript_name := chr19_tx[target_id == tid, transcript_name]]
            return(df1)
        }
    ))

    # limit to top 3 rows
    limit_tx <- df[transcript_id %in% tx_ids[1:3, target_id], unique(transcript_name)]

    # get condition and samples
    sample_cdn <- as.data.table(so_chr19$fits$full$design_matrix, keep.rownames = TRUE)

    gg <- (
        ggplot(
            data = df[transcript_name %in% limit_tx],
            mapping = aes(
                x = sample,
                fill = sample
            )
        )
        + geom_linerange(aes(
            ymin = min,
            ymax = max
        ))
        + geom_crossbar(aes(
            ymin = lower,
            y = mid,
            ymax = upper
        ))
        + scale_x_discrete(
            name = NULL,
            breaks = meta[, Sample_ID],
            labels = meta[, Label]
        )
        + scale_y_continuous(
            name = bquote("Estimated mRNA Count (" * log[e] * ")")
        )
        + scale_fill_manual(
            name = NULL,
            breaks = c(
                sample_cdn[conditionnovel == 1, rn],
                sample_cdn[conditionnovel == 0, rn]
            ),
            values = c(
                rep("#FF6347", sample_cdn[conditionnovel == 1, .N]),
                rep("#BDBDBD", sample_cdn[conditionnovel == 0, .N])
            )
        )
        + facet_grid(
            transcript_name ~ condition,
            scales = "free",
            labeller = as_labeller(labeller)
        )
        + guides(fill = FALSE)
        + theme_minimal()
        + theme(
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            axis.text.y = element_text(colour = "#000000"),
            panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_line(colour = "#000000"),
            panel.grid.minor.y = element_line(colour = "#000000")
        )
    )
    ggsave(
        file.path(PLOT_DIR, paste0(gene_name, ".dge.bootstraps.png")),
        gg,
        width = 16,
        height = 12,
        units = "cm"
    )
    ggsave(
        file.path(PLOT_DIR, paste0(gene_name, ".dge.bootstraps.pdf")),
        gg,
        width = 16,
        height = 12,
        units = "cm"
    )
}

gg_pval <- (
    ggplot(
        data = chr19_genes,
        mapping = aes(x = pval)
    )
    + geom_histogram(bins = 50)
    + theme_minimal()
)
ggsave(
    file.path(PLOT_DIR, "dge.chr19.png"),
    gg_pval,
    width = 16,
    height = 12,
    units = "cm"
)
ggsave(
    file.path(PLOT_DIR, "dge.chr19.pdf"),
    gg_pval,
    width = 16,
    height = 12,
    units = "cm"
)

gg_pval <- (
    ggplot(
        data = chrY_genes,
        mapping = aes(x = pval)
    )
    + geom_histogram(bins = 50)
    + theme_minimal()
)
ggsave(
    file.path(PLOT_DIR, "dge.chrY.png"),
    gg_pval,
    width = 16,
    height = 12,
    units = "cm"
)
ggsave(
    file.path(PLOT_DIR, "dge.chrY.pdf"),
    gg_pval,
    width = 16,
    height = 12,
    units = "cm"
)

gg_chr19_volcano <- (
    ggplot(
        data = chr19_tx,
        mapping = aes(
            x = b / log(2),
            y = -log10(pval),
            fill = paste(qval < 0.05, abs(b / log(2)) > 1, sign(b)),
            label = transcript_name
        )
    )
    + geom_point(shape = 21, alpha = 0.7)
    + geom_label_repel(
        data = subset(chr19_tx, qval < 0.05)
    )
    + scale_x_continuous(
        name = bquote(log[2] * "(FC)")
    )
    + scale_y_continuous(
        name = bquote(-log[10] * "(p-value)")
    )
    + scale_fill_manual(
        name = NULL,
        breaks = c(
            "FALSE FALSE 1",
            "FALSE FALSE -1",
            "FALSE TRUE 1",
            "FALSE TRUE -1",
            "TRUE TRUE 1", 
            "TRUE TRUE -1"
        ),
        labels = c(
            "N.S.",
            "N.S.",
            "N.S.",
            "N.S.",
            "q < 0.05, log2(FC) < -1",
            "q < 0.05, log2(FC) > 1"
        ),
        values = rep(
            c("#BDBDBD", "#FFA500", "#6495ED"),
            c(4, 1, 1)
        )
    )
    + guides(fill = FALSE)
    + theme_minimal()
)
ggsave(
    file.path(PLOT_DIR, "dge.chr19.png"),
    gg_chr19_volcano,
    width = 8,
    height = 6,
    units = "cm"
)
ggsave(
    file.path(PLOT_DIR, "dge.chr19.pdf"),
    gg_chr19_volcano,
    width = 8,
    height = 6,
    units = "cm"
)

gg_chrY_volcano <- (
    ggplot(
        data = chrY_tx,
        mapping = aes(
            x = b / log(2),
            y = -log10(pval),
            fill = paste(qval < 0.05, abs(b / log(2)) > 1, sign(b)),
            label = transcript_name
        )
    )
    + geom_point(shape = 21, alpha = 0.7)
    + geom_label_repel(
        data = subset(chrY_tx, qval < 0.05)
    )
        + scale_x_continuous(
        name = bquote(log[2] * "(FC)")
    )
    + scale_y_continuous(
        name = bquote(-log[10] * "(p-value)")
    )
    + scale_fill_manual(
        name = NULL,
        breaks = c(
            "FALSE FALSE 1",
            "FALSE FALSE -1",
            "FALSE TRUE 1",
            "FALSE TRUE -1",
            "TRUE TRUE 1", 
            "TRUE TRUE -1"
        ),
        labels = c(
            "N.S.",
            "N.S.",
            "N.S.",
            "N.S.",
            "q < 0.05, log2(FC) < -1",
            "q < 0.05, log2(FC) > 1"
        ),
        values = rep(
            c("#BDBDBD", "#FFA500", "#6495ED"),
            c(4, 1, 1)
        )
    )
    + guides(fill = FALSE)
    + theme_minimal()
)
ggsave(
    file.path(PLOT_DIR, "dge.chrY.png"),
    gg_chrY_volcano,
    width = 8,
    height = 6,
    units = "cm"
)
ggsave(
    file.path(PLOT_DIR, "dge.chrY.pdf"),
    gg_chrY_volcano,
    width = 8,
    height = 6,
    units = "cm"
)


# ==============================================================================
# Save Data
# ==============================================================================
loginfo("Saving data")

