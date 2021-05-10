# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))

# ==============================================================================
# Functions
# ==============================================================================
#' Save figures in multiple formats
#'
#' @param gg ggplot object
#' @param prefix Prefix for output file
#' @param ext Output extensions
#' @param dpi DPI resolution
savefig = function(gg, prefix, ext = c("png", "pdf"), width = 20, height = 12, dpi = 400) {
    for (e in ext) {
        ggsave(
            paste(prefix, e, sep = "."),
            gg,
            height = height,
            width = width,
            units = "cm",
            dpi = dpi
        )
    }
}

# ==============================================================================
# Data
# ==============================================================================
# load metadata
metadata <- fread(
    "../../Data/External/LowC_Samples_Data_Available.tsv",
    sep = "\t",
    header = TRUE,
)
# only keep included samples
metadata <- metadata[Include == "Yes"]
# SampleID column
metadata[, SampleID := paste0("Pca", get("Sample ID"))]

PLOT_DIR <- "Plots"

# ==============================================================================
# Analysis
# ==============================================================================
peaks = rbindlist(lapply(
    metadata[, SampleID],
    function(s) {
        print(s)
        dt = fread(
            paste0("../../Data/Processed/2019-05-03_PCa-H3K27ac-peaks/Peaks/", s, "_peaks.filtered.narrowPeak"),
            sep = "\t",
            select = c(1:3, 9),
            col.names = c("chr", "start", "end", "log10q")
        )
        dt[, SampleID := s]
        dt[, start := as.numeric(start)]
        dt[, end := as.numeric(end)]
        return(dt)
    }
))
peaks[, SampleID := as.factor(SampleID)]

peak_dists = rbindlist(lapply(
    metadata[, SampleID],
    function(s) {
        print(s)
        dt = fread(
            paste0("Closest/", s, ".closest-peaks.bed"),
            sep = "\t",
            select = 21,
            col.names = c("dist")
        )
        dt[, SampleID := s]
        return(dt)
    }
))
peak_dists[, SampleID := as.factor(SampleID)]

# ==============================================================================
# Plots
# ==============================================================================

# Number of peaks
# ------------------------------------------------------------------------------
n_peaks <- peaks[, .N, by = SampleID]
invisible(lapply(
    1:n_peaks[, .N],
    function(i) {
        t2e = metadata[SampleID == n_peaks[i, SampleID], get("T2E Status")]
        n_peaks[i, T2E_Status := t2e]
    }
))

gg_num_peaks <- (
    ggplot(data = n_peaks)
    + geom_col(aes(x = SampleID, y = N / 1000, fill = SampleID))
    + labs(x = NULL, y = expression("Number of Peaks (" * 10^3 * ")"))
    + scale_x_discrete(
        limits = metadata[order(get("Patient ID")), SampleID],
        labels = metadata[order(get("Patient ID")), get("Patient ID")]
    )
    + scale_fill_manual(
        limits = metadata[, SampleID],
        values = metadata[, Colour]
    )
    + guides(fill = FALSE)
    + theme_minimal()
    + theme(axis.text.x = element_text(angle = 90))
)
savefig(gg_num_peaks, file.path(PLOT_DIR, "peak-counts"))

# Peak size distribution
# ------------------------------------------------------------------------------
# reorder samples according to median peak size
median_peak_sizes = peaks[, median(end - start), by = SampleID]
sample_order = median_peak_sizes[order(V1), SampleID]
peaks[, SampleID := factor(SampleID, levels = sample_order, ordered  = TRUE)]

gg_peak_size <- (
    ggplot(data = peaks)
    + geom_violin(aes(x = SampleID, y = log10(end - start), fill = SampleID))
    + geom_boxplot(aes(x = SampleID, y = log10(end - start)), width = 0.2)
    + labs(x = NULL, y = expression(log[10] * "(Peak Size)"))
    + scale_x_discrete(
        limits = metadata[order(get("Patient ID")), SampleID],
        labels = metadata[order(get("Patient ID")), get("Patient ID")]
    )
    + scale_fill_manual(
        limits = metadata[, SampleID],
        values = metadata[, Colour]
    )
    + guides(fill = FALSE)
    + theme_minimal()
    + theme(axis.text.x = element_text(angle = 90))
)
savefig(gg_peak_size, file.path(PLOT_DIR, "peak-sizes"))

# Distance between
# ------------------------------------------------------------------------------
# reorder samples according to median peak size
median_dist_sizes = peak_dists[, median(dist), by = SampleID]
mean_dist_sizes = peak_dists[, mean(dist), by = SampleID]
sample_order = median_dist_sizes[order(V1), SampleID]
peak_dists[, SampleID := factor(SampleID, levels = sample_order, ordered  = TRUE)]

gg_peak_dist <- (
    ggplot(data = peak_dists)
    + geom_violin(aes(x = SampleID, y = log10(dist), fill = SampleID))
    + geom_boxplot(aes(x = SampleID, y = log10(dist)), width = 0.2)
    + labs(x = NULL, y = expression(log[10] * "(Distance (bp))"))
    + scale_x_discrete(
        limits = metadata[order(get("Patient ID")), SampleID],
        labels = metadata[order(get("Patient ID")), get("Patient ID")]
    )
    + scale_fill_manual(
        limits = metadata[, SampleID],
        values = metadata[, Colour]
    )
    + guides(fill = FALSE)
    + theme_minimal()
    + theme(axis.text.x = element_text(angle = 90))
)
savefig(gg_peak_dist, file.path(PLOT_DIR, "peak-dists"))


dist_centres = rbind(median_dist_sizes, mean_dist_sizes)
colnames(dist_centres) = c("SampleID", "Size")
dist_centres[, Type := rep(c("Median", "Mean"), each = metadata[, .N])]

gg_dist_centres <- (
    ggplot(data = dist_centres)
    + geom_col(aes(x = SampleID, y = Size, fill = SampleID))
    + labs(x = NULL, y = "Size (bp)", title = "H3K27ac centres of distance between peaks")
    + scale_x_discrete(
        limits = metadata[order(get("Patient ID")), SampleID],
        labels = metadata[order(get("Patient ID")), get("Patient ID")]
    )
    + scale_fill_manual(
        limits = metadata[, SampleID],
        values = metadata[, Colour]
    )
    + facet_wrap(. ~ Type, scales = "free_y")
    + guides(fill = FALSE)
    + theme_minimal()
    + theme(axis.text.x = element_text(angle = 90))
)
savefig(gg_peak_dist, file.path(PLOT_DIR, "peak-dists-centres"))
