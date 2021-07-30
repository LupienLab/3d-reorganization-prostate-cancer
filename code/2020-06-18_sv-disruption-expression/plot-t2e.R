# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("sleuth"))
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
source(file.path("..", "src", "savefig.R"))

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
        "..", "..", "Data", "External", "LowC_Samples_Data_Available.tsv"
    ),
    sep = "\t",
    header = TRUE
)
metadata <- metadata[Include == "Yes"]
metadata[, SampleID :=  paste0("PCa", get("Sample ID"))]
SAMPLES <- metadata$SampleID

so <- readRDS(file.path(RES_DIR, "sleuth-object.rds"))
so_transcripts <- fread(
    file.path(RES_DIR, "results.transcripts.tsv"),
    sep = "\t"
)
so_genes <- fread(
    file.path(RES_DIR, "results.genes.tsv"),
    sep = "\t"
)


# ==============================================================================
# Analysis
# ==============================================================================
# full table of estimated counts
full_table <- as.data.table(kallisto_table(so))

# erg-specific counts
erg <- rbindlist(lapply(
    so_transcripts[gene_name == "ERG", target_id],
    function(target) {
        dt <- as.data.table(get_bootstrap_summary(so, target, "est_counts"))
        dt[, target_id := target]
        return(dt)
    }
))

# ZNF516 counts
znf <- rbindlist(lapply(
    so_transcripts[gene_name == "ZNF516", target_id],
    function(target) {
        dt <- as.data.table(get_bootstrap_summary(so, target, "est_counts"))
        dt[, target_id := target]
        return(dt)
    }
))

erg_tpm <- merge(
    x = erg,
    y = full_table[, .(target_id, sample, tpm)],
    by = c("target_id", "sample")
)
colnames(erg_tpm)[colnames(erg_tpm) == "sample"] <- "SampleID"
colnames(erg_tpm)[colnames(erg_tpm) == "tpm"] <- "TPM"

znf_tpm <- merge(
    x = znf,
    y = full_table[, .(target_id, sample, tpm)],
    by = c("target_id", "sample")
)
colnames(znf_tpm)[colnames(znf_tpm) == "sample"] <- "SampleID"
colnames(znf_tpm)[colnames(znf_tpm) == "tpm"] <- "TPM"

# ==============================================================================
# Plots
# ==============================================================================
gg_erg <- (
    ggplot(data = erg_tpm)
    + geom_boxplot(
        aes(x = condition, y = TPM, fill = condition),
        alpha = 0.3,
        outlier.shape = NA
    )
    + geom_point(
        aes(x = condition, y = TPM, colour = condition),
        position = position_jitter(height = 0, width = 0.2)
    )
    + scale_x_discrete(
        breaks = c("No", "Yes"),
        labels = c("T2E-", "T2E+"),
        name = NULL
    )
    + guides(colour = FALSE, fill = FALSE)
    + facet_wrap(~ target_id, drop = TRUE, scale = "free")
    + theme_minimal()
    + theme(
        axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 1)
    )
)
savefig(
    gg_erg,
    file.path(PLOT_DIR, "ERG-expression"),
    height = 20,
    width = 20
)

gg_erg_sum <- (
    ggplot(data = erg_tpm)
    + geom_boxplot(
        aes(x = condition, y = sum(TPM), fill = condition),
        alpha = 0.3,
        outlier.shape = NA
    )
    + geom_point(
        aes(x = condition, y = sum(TPM), colour = condition),
        position = position_jitter(height = 0, width = 0.2)
    )
    + scale_x_discrete(
        breaks = c("No", "Yes"),
        labels = c("T2E-", "T2E+"),
        name = NULL
    )
    + guides(colour = FALSE, fill = FALSE)
    + theme_minimal()
    + theme(
        axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 1)
    )
)
savefig(
    gg_erg_sum,
    file.path(PLOT_DIR, "ERG-expression.gene-level")
)

gg_znf_sum <- (
    ggplot(data = znf_tpm[, .(TPM = sum(TPM)), by = c("condition", "SampleID")])
    + geom_boxplot(
        aes(x = condition, y = TPM, fill = condition),
        alpha = 0.3,
        outlier.shape = NA
    )
    + geom_point(
        aes(x = condition, y = TPM, colour = condition),
        position = position_jitter(height = 0, width = 0.2)
    )
    + geom_path(
        data = data.table(
            x = rep(c("No", "Yes"), each = 2),
            y = c(10, 25, 25, 24.5),
            group = 1
        ),
        mapping = aes(x = x, y = y, group = group)
    )
    + geom_text(
        data = data.table(
            x = 1.5,
            y = 25,
            label = paste0("p = ", so_genes[gene_name == "ZNF516", signif(pval, digits = 3)])
        ),
        mapping = aes(x = x, y = y, label = label),
        vjust = -0.8
    )
    + scale_x_discrete(
        breaks = c("No", "Yes"),
        labels = c("T2E-", "T2E+"),
        name = NULL
    )
    + scale_colour_manual(
        breaks = c("No", "Yes"),
        labels = c("T2E-", "T2E+"),
        values = c("#418B3D", "#3215C1"),
        name = NULL
    )
    + guides(colour = FALSE, fill = FALSE)
    + theme_minimal()
    + theme(
        axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 1)
    )
)
savefig(
    gg_znf_sum,
    file.path(PLOT_DIR, "ZNF516-expression.gene-level"),
    width = 8
)
