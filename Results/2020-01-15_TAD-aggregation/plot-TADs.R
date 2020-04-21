# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))

MAX_WINDOW <- 24
MIN_WINDOW <- 3
MAX_PERSISTENCE <- MAX_WINDOW - MIN_WINDOW + 1
PLOT_DIR <- "Plots"

# ==============================================================================
# Functions
# ==============================================================================
split_comma_col <- function(v, f=identity) {
    # split into list
    splitv <- lapply(v, function(x) {strsplit(x, "\\|")[[1]]})
    # remove various non-informative characters (spaces, braces)
    splitv <- lapply(splitv, function(x) {gsub("[][ ]", "", x)})
    return(lapply(splitv, f))
}

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
metadata <- fread(file.path("..", "..", "Data", "External", "LowC_Samples_Data_Available.tsv"))
metadata <- metadata[Include == "Yes"]
metadata[, SampleID := paste0("PCa", get("Sample ID"))]
SAMPLES <- metadata[, SampleID]

# load aggregated boundary calls from each sample
boundaries = rbindlist(lapply(
    SAMPLES,
    function(s) {
        dt = fread(
            file.path("resolved-TADS", paste0(s, ".40000bp.aggregated-boundaries.tsv")),
            sep = "\t",
            header = TRUE
        )
        dt[, SampleID := s]
        return(dt)
    }
))
boundaries$w <- split_comma_col(boundaries$w, as.numeric)

# load TADs
tads <- rbindlist(lapply(
    SAMPLES,
    function(s) {
        dt <- fread(paste0("resolved-TADs/", s, ".40000bp.aggregated-domains.bedGraph"),
            sep = "\t",
            header = FALSE,
            col.names = c("chr", "start", "end", "lower_persistence", "upper_persistence", "w")
        )
        dt[, SampleID := s]
    }
))
tads[, width := as.numeric(end - start)]

# only keep TADs called at w <= MAX_WINDOW
tads <- tads[w <= MAX_WINDOW]
tads[, lower_persistence := pmin(MAX_PERSISTENCE, lower_persistence)]
tads[, upper_persistence := pmin(MAX_PERSISTENCE, upper_persistence)]
boundaries[, Persistence := pmin(MAX_PERSISTENCE, Order)]

# add Patient IDs to tads and boundaries, not just SampleIDs
tads <- merge(tads, metadata[, .SD, .SDcols = c("SampleID", "Patient ID")])
colnames(tads)[which(colnames(tads) == "Patient ID")] <- "Patient_ID"
boundaries <- merge(boundaries, metadata[, .SD, .SDcols = c("SampleID", "Patient ID")])
colnames(boundaries)[which(colnames(boundaries) == "Patient ID")] <- "Patient_ID"

# load BPscore calculations
bpscore <- fread("Statistics/tad-distances.tsv", sep = "\t", header = TRUE)
window_diffs <- fread("Statistics/tad-similarity-deltas.tsv", sep = "\t", header = TRUE)

# load CTCF distances
ctcf_pairs = fread("CTCF/TAD-boundary.LNCaP-CTCF-peaks.distances.tsv", sep = "\t", header = TRUE)

# ==============================================================================
# Analysis
# ==============================================================================
boundary_counts = boundaries[, .N, by = c("SampleID", "Patient_ID")]
boundary_counts_persistence = boundaries[, .N, by = c("SampleID", "Patient_ID", "Persistence")]

# calculate coefficient of variation across TAD sizes to see where samples vary
tads_mean_size <- tads[, mean(width), by = c("SampleID", "Patient_ID", "w")]
tads_cov <- tads_mean_size[, sd(V1) / mean(V1), by = "w"]
tads_median_size <- tads[, median(width), by = c("SampleID", "Patient_ID", "w")]
tads_madm <- tads_median_size[, median(abs(V1 - median(V1))), by = "w"]

# calculate CTCF peak proximity to TAD boundaries
ctcf_pairs[, Distance_Bin := sign(Distance) * floor(abs(Distance) / 1e3)]
ctcf_pairs_binned <- ctcf_pairs[,
    .N,
    by = c("SampleID", "Distance_Bin", "chr_bound", "start_bound", "end_bound")
][, .(Mean_N = mean(N)), keyby = c("SampleID", "Distance_Bin")]

# ==============================================================================
# Plots
# ==============================================================================
# Boundaries
# --------------------------------------
# plot number of resolved boundaries
gg_boundaries <- (
    ggplot(data = boundary_counts)
    + geom_col(aes(x = SampleID, y = N, fill = SampleID))
    + scale_x_discrete(
        limits = metadata[, SampleID],
        labels = metadata[, get("Patient ID")]
    )
    + scale_fill_manual(
        limits = metadata[, SampleID],
        values = metadata[, Colour]
    )
    + labs(x = NULL, y = "Number of unique boundaries")
    + guides(fill = FALSE)
    + theme_minimal()
    + theme(
        axis.text.x = element_text(angle = 90)
    )
)
savefig(gg_boundaries, file.path(PLOT_DIR, "boundary-counts"))

# plot number of resolved boundaries by order
gg_bounds_persistence <- (
    ggplot(data = boundary_counts_persistence)
    + geom_col(aes(x = Persistence, y = N, fill = SampleID, group = SampleID), position = "dodge")
    + labs(x = "Boundary Persistence", y = "Number of unique boundaries")
    + scale_x_discrete(
        breaks = c(1, 6, 11, 16, MAX_PERSISTENCE),
        limits = c(1, 6, 11, 16, MAX_PERSISTENCE),
        labels = c(1, 6, 11, 16, MAX_PERSISTENCE)
    )
    + scale_fill_manual(
        limits = metadata[, SampleID],
        labels = metadata[, get("Patient ID")],
        values = metadata[, Colour],
        name = "Patient"
    )
    + theme_minimal()
    + theme(
        axis.text.x = element_text(angle = 90)
    )
)
savefig(gg_bounds_persistence, file.path(PLOT_DIR, "boundary-counts.by-persistence"))

# CTCF binding site proximity to boundaries
gg_bounds_ctcf <- (
    ggplot(data = ctcf_pairs)
    + stat_density(
        aes(x = Distance / 1e3, y = ..count.., colour = SampleID),
        bw = 1,
        geom = "path",
        position = position_identity()
    )
    + labs(x = "Distance from TAD boundary (kbp)", y = "LNCaP CTCF Peaks")
    + scale_x_continuous(
        limits = c(-150, 150),
        breaks = seq(-150, 150, 50),
        labels = seq(-150, 150, 50),
    )
    + scale_colour_manual(
        limits = metadata[, SampleID],
        labels = metadata[, get("Patient ID")],
        values = metadata[, Colour],
        name = "Patient"
    )
    + theme_minimal()
)
# gg_bounds_ctcf <- (
#     ggplot(data = ctcf_pairs_binned)
#     + geom_path(aes(x =  Distance_Bin, y = Mean_N, colour = SampleID, group = SampleID))
#     + labs(x = "Distance from TAD boundary (kbp)", y = "LNCaP CTCF Peaks / 5 kbp")
#     + scale_x_continuous(
#         limits = c(-150, 150),
#         breaks = seq(-150, 150, 50),
#         labels = seq(-150, 150, 50),
#     )
#     + scale_colour_manual(
#         limits = metadata[, SampleID],
#         labels = metadata[, get("Patient ID")],
#         values = metadata[, Colour],
#         name = "Patient"
#     )
#     + theme_minimal()
# )
savefig(gg_bounds_ctcf, file.path(PLOT_DIR, "boundary-counts.ctcf-proximity"))

# TADs
# --------------------------------------
# TAD counts per sample faceted by window size
gg_tad_counts_window <- (
    ggplot(data = tads[, .N, by = c("Patient_ID", "w")])
    + geom_col(aes(x = Patient_ID, y = N, fill = Patient_ID))
    + labs(x = NULL, y = "Number of TADs")
    + scale_fill_manual(
        limits = metadata[, get("Patient ID")],
        labels = metadata[, get("Patient ID")],
        values = metadata[, Colour]
    )
    + guides(fill = FALSE)
    + facet_wrap(~ w, nrow = 6)
    + theme_minimal()
    + theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "bottom"
    )
)
savefig(gg_tad_counts_window, file.path(PLOT_DIR, "tad-counts.by-window"), height = 20)

# TAD counts per sample faceted by window size
gg_tad_counts_sample <- (
    ggplot(data = tads[, .N, by = c("Patient_ID", "w")])
    + geom_col(aes(x = w, y = N, fill = w))
    + labs(x = "Window size", y = "Number of TADs")
    + scale_fill_viridis_c()
    + guides(fill = FALSE)
    + facet_wrap(~ Patient_ID, nrow = 3)
    + theme_minimal()
    + theme(
        legend.position = "bottom"
    )
)
savefig(gg_tad_counts_sample, file.path(PLOT_DIR, "tad-counts.by-sample"), height = 20)

# TAD counts per sample at each window size
gg_tad_counts <- (
    ggplot(data = tads[, .N, by = c("Patient_ID", "w")])
    + geom_path(aes(x = w, y = N, colour = Patient_ID, group = Patient_ID))
    + geom_point(aes(x = w, y = N, colour = Patient_ID))
    + labs(x = "Window size", y = "Number of TADs")
    + scale_x_discrete(
        breaks = seq(MIN_WINDOW, MAX_WINDOW, by = 3),
        limits = seq(MIN_WINDOW, MAX_WINDOW, by = 3),
        labels = seq(MIN_WINDOW, MAX_WINDOW, by = 3)
    )
    + scale_colour_manual(
        limits = metadata[, get("Patient ID")],
        values = metadata[, Colour],
        name = "Patient"
    )
    + theme_minimal()
)
savefig(gg_tad_counts, file.path(PLOT_DIR, "tad-counts"))

# TAD size distribution
gg_tad_size <- (
    ggplot(data = tads)
    + stat_density(
        aes(x = width / 1e6, y = 0.1 * ..count.., colour = SampleID),
        alpha = 0.8,
        geom = "path",
        position = position_identity()
    )
    + labs(x = "TAD Size (Mbp)", y = "Number of TADs")
    + scale_x_continuous(
        limits = c(0, 5)
    )
    + scale_colour_manual(
        limits = metadata[, SampleID],
        labels = metadata[, get("Patient ID")],
        values = metadata[, Colour],
        name = "Patient"
    )
    + facet_wrap(~ w, scales = "free_y", nrow = 5)
    + theme_minimal()
    + theme(
        legend.position = "bottom"
    )
)
savefig(gg_tad_size, file.path(PLOT_DIR, "tad-size.distribution"), height = 20)

# same as above but only focussing on a few window sizes
gg_tad_size_reduced <- (
    ggplot(data = tads[w %in% c(MIN_WINDOW, floor((MAX_WINDOW + MIN_WINDOW) / 2), MAX_WINDOW)])
    + stat_density(
        aes(x = width / 1e6, y = 0.1 * ..count.., colour = SampleID),
        alpha = 0.8,
        geom = "path",
        position = position_identity()
    )
    + labs(x = "TAD Size (Mbp)", y = "Number of TADs")
    + scale_x_continuous(
        limits = c(0, 5)
    )
    + scale_colour_manual(
        limits = metadata[, SampleID],
        labels = metadata[, get("Patient ID")],
        values = metadata[, Colour],
        name = "Patient"
    )
    + facet_wrap(~ w, ncol = 3, scales = "free_y", labeller = label_both)
    + theme_minimal()
    + theme()
)
savefig(gg_tad_size_reduced, file.path(PLOT_DIR, "tad-size.distribution.reduced"))

# Variance between samples
# --------------------------------------
# coefficient of variation for BPscore
gg_cov <- (
    ggplot(data = tads_cov)
    + geom_col(aes(x = w, y = V1, fill = w))
    + labs(x = "Window size", y = "Coefficient of Variation")
    + scale_x_discrete(
        breaks = seq(MIN_WINDOW, MAX_WINDOW, by = 3),
        limits = seq(MIN_WINDOW, MAX_WINDOW, by = 3),
        labels = seq(MIN_WINDOW, MAX_WINDOW, by = 3)
    )
    + scale_fill_viridis_c()
    + guides(fill = FALSE)
    + theme_minimal()
)
savefig(gg_cov, file.path(PLOT_DIR, "tad-size.coeff-var"))

# plot distances between TADs of different samples across window sizes
gg_bpscore = (
    ggplot(data = bpscore[w <= MAX_WINDOW])
    + geom_point(aes(x = w, y = 1 - dist, colour = w), position = position_jitter(width = 0.2, height = 0))
    + geom_boxplot(aes(x = w, y = 1 - dist, group = w), alpha = 0.5, outlier.shape = NA, width = 0.5)
    + labs(x = "Window size", y = "TAD similarity (1 - BPscore)")
    + scale_x_discrete(
        breaks = seq(MIN_WINDOW, MAX_WINDOW, by = 3),
        limits = seq(MIN_WINDOW, MAX_WINDOW, by = 3),
        labels = seq(MIN_WINDOW, MAX_WINDOW, by = 3)
    )
    + scale_colour_viridis_c()
    + scale_fill_viridis_c()
    + guides(colour = FALSE)
    + theme_minimal()
)
savefig(gg_bpscore, file.path(PLOT_DIR, "bp-score"))

# plot the change in similarities over window sizes
gg_sim_window = (
    ggplot(data = window_diffs[w <= MAX_WINDOW])
    + geom_path(aes(x = w, y = 1 - diff, group = SampleID, colour = SampleID))
    + labs(x = "Window size", y = expression("1 - " * delta["w, w-1"]))
    + scale_x_discrete(
        breaks = seq(MIN_WINDOW, MAX_WINDOW, by = 3),
        limits = seq(MIN_WINDOW, MAX_WINDOW, by = 3),
        labels = seq(MIN_WINDOW, MAX_WINDOW, by = 3)
    )
    + scale_colour_manual(
        limits = metadata[, SampleID],
        labels = metadata[, get("Patient ID")],
        values = metadata[, Colour],
        name = "Patient"
    )
    + theme_minimal()
)
savefig(gg_sim_window, file.path(PLOT_DIR, "bp-score.window-similarity"))


# plot the change in similarities over window sizes
gg_sim_window_delta = (
    ggplot(data = window_diffs[w >= 5 & w <= MAX_WINDOW])
    + geom_path(aes(x = w, y = abs_delta, group = SampleID, colour = SampleID))
    + labs(x = "Window size", y = expression("|" * delta["w, w-1"] - delta["w-1, w-2"] * "|"))
    + scale_x_discrete(
        breaks = seq(MIN_WINDOW, MAX_WINDOW, by = 3),
        limits = seq(MIN_WINDOW, MAX_WINDOW, by = 3),
        labels = seq(MIN_WINDOW, MAX_WINDOW, by = 3)
    )
    + scale_colour_manual(
        limits = metadata[, SampleID],
        labels = metadata[, get("Patient ID")],
        values = metadata[, Colour],
        name = "Patient"
    )
    + theme_minimal()
)
savefig(gg_sim_window_delta, file.path(PLOT_DIR, "bp-score.window-similarity.delta"))
