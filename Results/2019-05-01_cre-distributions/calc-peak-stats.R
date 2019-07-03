# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))

# ==============================================================================
# Data
# ==============================================================================
metadata = fread(
    "../../Data/External/LowC_Samples_Data_Available.tsv",
    sep = "\t",
    header = TRUE,
    select = c(1, 3),
    col.names = c("SampleID", "T2E_Status")
)
metadata[, SampleID := paste0("Pca", SampleID)]

# ==============================================================================
# Analysis
# ==============================================================================
peaks = rbindlist(lapply(
    metadata[, SampleID],
    function(s) {
        print(s)
        dt = fread(
            paste0("../../Data/Processed/2019-05-03_PCa-H3K27ac-peaks/", s, "_peaks.filtered.narrowPeak"),
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
n_peaks = peaks[, .N, by = SampleID]
sample_order = n_peaks[order(N), as.character(SampleID)]
n_peaks[, SampleID := factor(SampleID, levels = sample_order, ordered  = TRUE)]

gg = (
    ggplot(data = n_peaks)
    + geom_col(aes(x = SampleID, y = N / 1000, fill = SampleID))
    + labs(x = "Patient Sample", y = expression("Number of Peaks (" * 10^3 * ")"), title = "H3K27ac peak counts")
    + guides(fill = FALSE)
    + theme(axis.text.x = element_text(angle = 90))
)
ggsave(
    "Stats/peak-counts.png",
    height = 12,
    width = 20,
    units = "cm"
)

# Peak size distribution
# ------------------------------------------------------------------------------
# reorder samples according to median peak size
median_peak_sizes = peaks[, median(end - start), by = SampleID]
sample_order = median_peak_sizes[order(V1), SampleID]
peaks[, SampleID := factor(SampleID, levels = sample_order, ordered  = TRUE)]

gg = (
    ggplot(data = peaks)
    + geom_violin(aes(x = SampleID, y = log10(end - start), fill = SampleID))
    + geom_boxplot(aes(x = SampleID, y = log10(end - start)), width = 0.2)
    + labs(x = "Patient Sample", y = "log10(Peak Size)", title = "H3K27ac peak sizes")
    + guides(fill = FALSE)
    + theme(axis.text.x = element_text(angle = 90))
)
ggsave(
    "Stats/peak-sizes.png",
    height = 12,
    width = 20,
    units = "cm"
)

# Distance between
# ------------------------------------------------------------------------------
# reorder samples according to median peak size
median_dist_sizes = peak_dists[, median(dist), by = SampleID]
mean_dist_sizes = peak_dists[, mean(dist), by = SampleID]
sample_order = median_dist_sizes[order(V1), SampleID]
peak_dists[, SampleID := factor(SampleID, levels = sample_order, ordered  = TRUE)]

gg = (
    ggplot(data = peak_dists)
    + geom_violin(aes(x = SampleID, y = log10(dist), fill = SampleID))
    + geom_boxplot(aes(x = SampleID, y = log10(dist)), width = 0.2)
    + labs(x = "Patient Sample", y = "log10(Distance)", title = "H3K27ac distance between peaks")
    + guides(fill = FALSE)
    + theme(axis.text.x = element_text(angle = 90))
)
ggsave(
    "Stats/peak-dists.png",
    height = 12,
    width = 20,
    units = "cm"
)

dist_centres = rbind(median_dist_sizes, mean_dist_sizes)
colnames(dist_centres) = c("SampleID", "Size")
dist_centres[, Type := rep(c("Median", "Mean"), each = metadata[, .N])]

gg = (
    ggplot(data = dist_centres)
    + geom_col(aes(x = SampleID, y = Size, fill = SampleID))
    + labs(x = "Patient Sample", y = "Size (bp)", title = "H3K27ac centres of distance between peaks")
    + facet_grid(. ~ Type, scales = "free")
    + guides(fill = FALSE)
    + theme(axis.text.x = element_text(angle = 90))
)
ggsave(
    "Stats/peak-dists-centres.png",
    height = 12,
    width = 20,
    units = "cm"
)
