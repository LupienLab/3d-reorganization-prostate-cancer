# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("sleuth"))
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
source("../2020-02-19_chromoplexy/plotting-helper.R")


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

so <- readRDS("sleuth-object.rds")
so_transcripts <- fread("t2e-comparison.transcripts.tsv", sep = "\t")


# ==============================================================================
# Analysis
# ==============================================================================
erg <- rbindlist(lapply(
    so_transcripts[gene_name == "ERG", target_id],
    function(target) {
        dt <- as.data.table(get_bootstrap_summary(so, target, "est_counts"))
        dt[, target_id := target]
        return(dt)
    }
))

full_table <- as.data.table(kallisto_table(so))

erg_tpm <- merge(
    x = erg,
    y = full_table[, .(target_id, sample, tpm)],
    by = c("target_id", "sample")
)
colnames(erg_tpm)[colnames(erg_tpm) == "sample"] <- "SampleID"
colnames(erg_tpm)[colnames(erg_tpm) == "tpm"] <- "TPM"

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
savefig(gg_erg, "Plots/ERG-expression", height = 20, width = 20)
