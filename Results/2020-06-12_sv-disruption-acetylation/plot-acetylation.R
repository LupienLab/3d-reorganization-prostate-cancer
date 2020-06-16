# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("pheatmap"))
suppressMessages(library("RColorBrewer"))
suppressMessages(library("ggplot2"))
source("../2020-02-19_chromoplexy/plotting-helper.R")

PLOT_DIR <- "Plots"

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

# load TAD acetylation hypothesis test values
acetyl <- fread(
    "sv-disruption-tests.acetylation.tsv",
    sep = "\t",
    header = TRUE
)

size_factors <- fread(
    "sv-disruption-tests.acetylation.size-factors.tsv",
    sep = "\t",
    header = TRUE
)

test_group_acetyl <- rbindlist(lapply(
    acetyl[, test_ID],
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
    filename = file.path(PLOT_DIR, "H3K27ac-correlation-tad-induced.png")
)

gg_qq <- (
    ggplot(data = acetyl[order(p)])
    + geom_point(aes(x = -log10(ppoints(acetyl[, .N])), y = -log10(p)))
    + geom_abline(aes(intercept = 0, slope = 1), linetype = "dashed")
    + scale_x_continuous(
        limits = c(0, 3)
    )
    + scale_y_continuous(
        limits = c(0, 15)
    )
    + labs(
        x = expression(-log[10] * "(Uniform quantile)"),
        y = expression(-log[10] * "(Observed quantile)")
    )
    + theme_minimal()
)
savefig(gg_qq, file.path(PLOT_DIR, "acetylation.qq-plot"))

gg_pval_hist <- (
    ggplot(data = acetyl)
    + geom_histogram(aes(x = p))
    + labs(x = "p-value", y = "Frequency")
    + theme_minimal()
)
savefig(gg_pval_hist, file.path(PLOT_DIR, "acetylation.p-values"))


gg_volcano <- (
    ggplot(data = acetyl)
    + geom_point(aes(x = z, y = -log10(padj)))
    + geom_hline(aes(yintercept = -log10(0.05)), linetype = "dashed")
    + labs(x = "Stouffer's Z", y = expression(-log[10] * " FDR"))
    + theme_minimal()
)
savefig(gg_volcano, file.path(PLOT_DIR, "acetylation.volcano"))

gg_all_tests_pval_hist <- (
    ggplot(data = test_group_acetyl)
    + geom_density(aes(x = p, colour = test_ID, group = test_ID))
    + labs(x = "p-value", y = "Frequency")
    + theme_minimal()
)
savefig(gg_all_tests_pval_hist, file.path(PLOT_DIR, "acetylation.p-values.all-tests"))

for (tid in acetyl[, test_ID]) {
    print(tid)
    gg_test_pval_hist <- (
        ggplot(data = test_group_acetyl[test_ID == tid, .SD])
        + geom_density(aes(x = p))
        + labs(x = "p-value", y = "Density")
        + theme_minimal()
    )
    savefig(gg_test_pval_hist, file.path(PLOT_DIR, "Tests", paste0("acetylation.p-values.test_", tid)))
}
