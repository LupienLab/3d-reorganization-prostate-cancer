# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("pheatmap"))
suppressMessages(library("RColorBrewer"))
suppressMessages(library("ggplot2"))
# ==============================================================================
# Data
# ==============================================================================
# load sample metadata
metadata <- fread(
    "../../Data/External/LowC_Samples_Data_Available.tsv",
    sep = "\t",
    header = TRUE
)
SAMPLES <- paste0("PCa", metadata[, get("Sample ID")])
metadata[, ChIP_file := paste0("../../Data/Processed/2019-05-03_PCa-H3K27ac-peaks/BAMs/Pca", get("Sample ID"), "_H3K27ac.sorted.dedup.bam")]
metadata[, Ctrl_file := paste0("../../Data/Processed/2019-05-03_PCa-H3K27ac-peaks/BAMs/Pca", get("Sample ID"), "_input.sorted.dedup.bam")]

# load raw counts of H3K27ac pull down and input, then take differences and get total ChIP library sizes
library_sizes <- sapply(SAMPLES, function(s) 0)
for (s in SAMPLES) {
    dt_chip <- fread(
        file.path("Acetylation", paste0(s, "_H3K27ac.induced-region-counts.bed")),
        sep = "\t",
        header = FALSE,
        col.names = c("chr", "start", "end", "count", "supported", "width", "frac_supported")
    )
    dt_input <- fread(
        file.path("Acetylation", paste0(s, "_input.induced-region-counts.bed")),
        sep = "\t",
        header = FALSE,
        col.names = c("chr", "start", "end", "count", "supported", "width", "frac_supported")
    )
    # calculate the linear scaling based on library size
    library_sizes[s] <- min(sample_library_size["chip"], sample_library_size["input"])
}

# load acetylation correlations
sample_corrs <- as.matrix(fread("acetylation-correlation.tsv", sep = "\t", header = TRUE))
rownames(sample_corrs) <- SAMPLES

# load TAD acetylation hypothesis test values
tads <- fread(
    "sv-disruption-tests.acetylation.tsv",
    sep = "\t",
    header = TRUE
)

# ==============================================================================
# Plots
# ==============================================================================
# heatmap of acetylation across the genome for all samples
ann_cols <- data.frame(
    T2E = metadata[, get("T2E Status")],
    LibrarySize = library_sizes
)
rownames(ann_cols) <- SAMPLES
pheatmap(
    mat = sample_corrs,
    color = colorRampPalette(brewer.pal(n = 7, name = "YlOrRd"))(100),
    #breaks = seq(0.8, 1, 0.01),
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    clustering_method = "ward.D2",
    annotation_col = ann_cols,
    legend = TRUE,
    filename = "Plots/H3K27ac-correlation-tad-induced.png"
)

gg = (
    ggplot(data = tads[order(p)])
    + geom_point(aes(x = -log10(ppoints(tads[, .N])), y = -log10(p)))
    + geom_abline(aes(intercept = 0, slope = 1), linetype = "dashed")
    + labs(
        x = expression(-log[10] * "(Uniform quantile)"),
        y = expression(-log[10] * "(Observed quantile)")
    )
    + theme_minimal()
)
ggsave(
    "Plots/sv-disruption.acetylation.qq-plot.png",
    height = 12,
    width = 20,
    units = "cm"
)

gg = (
    ggplot(data = tads)
    + geom_histogram(aes(x = p))
    + labs(x = "p-value", y = "Frequency")
    + theme_minimal()
)
ggsave(
    "Plots/sv-disruption.acetylation.p-values.png",
    height = 12,
    width = 20,
    units = "cm"
)


gg = (
    ggplot(data = tads)
    + geom_point(aes(x = z, y = -log10(padj)))
    + geom_hline(aes(yintercept = -log10(0.05)), linetype = "dashed")
    + labs(x = "Stouffer's Z", y = expression(-log[10] * " FDR"))
    + theme_minimal()
)
ggsave(
    "Plots/sv-disruption.acetylation.volcano.png",
    height = 12,
    width = 20,
    units = "cm"
)
