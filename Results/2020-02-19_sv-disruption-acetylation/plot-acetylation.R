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
metadata <- metadata[Include == "Yes"]
metadata[, SampleID := paste0("PCa", get("Sample ID"))]
metadata[, ChIP_file := paste0("../../Data/Processed/2019-05-03_PCa-H3K27ac-peaks/BAMs/Pca", get("Sample ID"), "_H3K27ac.sorted.dedup.bam")]
metadata[, Ctrl_file := paste0("../../Data/Processed/2019-05-03_PCa-H3K27ac-peaks/BAMs/Pca", get("Sample ID"), "_input.sorted.dedup.bam")]
SAMPLES <- metadata[, SampleID]

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
}
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
    library_sizes[s] <- min(dt_chip[, sum(count)], dt_input[, sum(count)])
}
 
# load acetylation correlations
sample_corrs <- as.matrix(fread("acetylation-correlation.tsv", sep = "\t", header = TRUE))
rownames(sample_corrs) <- SAMPLES

# load TAD acetylation hypothesis test values
acetyl <- fread(
    "Graphs/sv-disruption-tests.acetylation.tsv",
    sep = "\t",
    header = TRUE
)

tests <- fread(
    "Graphs/sv-disruption-tests.tsv",
    sep = "\t",
    header = TRUE
)

test_group_acetyl <- rbindlist(lapply(
    tests[, test_ID],
    function(tid) {
        dt <- fread(
            paste0("Acetylation/Tests/test_", tid, ".results.tsv"),
            sep = "\t",
            header = TRUE
        )
        dt[, test_ID := tid]
        return(dt)
    }
))
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

gg_qq = (
    ggplot(data = acetyl[order(p)])
    + geom_point(aes(x = -log10(ppoints(acetyl[, .N])), y = -log10(p)))
    + geom_abline(aes(intercept = 0, slope = 1), linetype = "dashed")
    + labs(
        x = expression(-log[10] * "(Uniform quantile)"),
        y = expression(-log[10] * "(Observed quantile)")
    )
    + theme_minimal()
)
ggsave(
    "Plots/sv-disruption/acetylation.qq-plot.png",
    gg_qq,
    height = 12,
    width = 20,
    units = "cm"
)
ggsave(
    "Plots/sv-disruption/acetylation.qq-plot.pdf",
    gg_qq,
    height = 12,
    width = 20,
    units = "cm"
)

gg_pval_hist = (
    ggplot(data = acetyl)
    + geom_histogram(aes(x = p))
    + labs(x = "p-value", y = "Frequency")
    + theme_minimal()
)
ggsave(
    "Plots/sv-disruption/acetylation.p-values.png",
    gg_pval_hist,
    height = 12,
    width = 20,
    units = "cm"
)
ggsave(
    "Plots/sv-disruption/acetylation.p-values.pdf",
    gg_pval_hist,
    height = 12,
    width = 20,
    units = "cm"
)


gg_volcano = (
    ggplot(data = acetyl)
    + geom_point(aes(x = z, y = -log10(padj)))
    + geom_hline(aes(yintercept = -log10(0.05)), linetype = "dashed")
    + labs(x = "Stouffer's Z", y = expression(-log[10] * " FDR"))
    + theme_minimal()
)
ggsave(
    "Plots/sv-disruption/acetylation.volcano.png",
    gg_volcano,
    height = 12,
    width = 20,
    units = "cm"
)
ggsave(
    "Plots/sv-disruption/acetylation.volcano.pdf",
    gg_volcano,
    height = 12,
    width = 20,
    units = "cm"
)

gg_all_tests_pval_hist <- (
    ggplot(data = test_group_acetyl)
    + geom_density(aes(x = p, colour = test_ID, group = test_ID))
    + labs(x = "p-value", y = "Frequency")
    + theme_minimal()
)
ggsave(
    "Plots/sv-disruption/acetylation.p-values.all-tests.png",
    gg_all_tests_pval_hist,
    height = 12,
    width = 20,
    units = "cm"
)

for (tid in tests[, test_ID]) {
    print(tid)
    gg_test_pval_hist <- (
        ggplot(data = test_group_acetyl[test_ID == tid, .SD])
        + geom_density(aes(x = p))
        + labs(x = "p-value", y = "Density")
        + theme_minimal()
    )
    ggsave(
        paste0("Acetylation/Tests/p-values.", tid, ".png"),
        gg_test_pval_hist,
        height = 12,
        width = 20,
        units = "cm"
    )
}
