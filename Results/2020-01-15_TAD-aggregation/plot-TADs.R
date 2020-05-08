# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
suppressMessages(library("scales"))
suppressMessages(library("MASS"))

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

# ECDFs number of TADs of a given width
tad_size_ecdf <- rbindlist(lapply(
    SAMPLES,
    function(s) {
        rbindlist(lapply(
            seq(MIN_WINDOW, MAX_WINDOW, 1),
            function(window) {
                f <- ecdf(tads[SampleID == s & w == window, width])
                return(data.table(
                    SampleID = s,
                    w = window,
                    x = seq(0, 5e6, 4e4),
                    y = f(seq(0, 5e6, 4e4))
                ))
            }
        ))
    }
))

tad_size_ecdf_est <- tad_size_ecdf[, .(Mean = mean(y), SD = sd(y)), by = c("w", "x")]
tad_size_ecdf_est[, Lower := Mean - SD]
tad_size_ecdf_est[, Upper := Mean + SD]

ctcf_fc <- ctcf_pairs[Bin_Mid == 0, .(Peak=Freq), keyby = "SampleID"]
ctcf_fc$Background <- ctcf_pairs[abs(Bin_Mid) > 100000, mean(Freq), keyby = "SampleID"]$V1

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
    + geom_path(aes(x = Bin_Mid / 1e3, y = Freq, colour = SampleID, group = SampleID))
    + labs(x = "Distance from TAD boundary (kbp)", y = "Average # LNCaP CTCF Peaks / 5 kbp")
    + scale_x_continuous(
        breaks = seq(-150, 150, 50),
        labels = seq(-150, 150, 50),
    )
    + scale_y_continuous(
        limits = c(0, 0.025)
    )
    + scale_colour_manual(
        limits = metadata[, SampleID],
        labels = metadata[, get("Patient ID")],
        values = metadata[, Colour],
        name = "Patient"
    )
    + coord_cartesian(xlim = c(-150, 150))
    + theme_minimal()
)
savefig(gg_bounds_ctcf, file.path(PLOT_DIR, "boundary-counts.ctcf-proximity"))

gg_bounds_ctcf_fc <- (
    ggplot(data = ctcf_fc)
    + geom_col(aes(x = SampleID, y = Peak / Background, fill = SampleID))
    + labs(x = NULL, y = "Fold change (peak vs background)")
    + scale_x_discrete(
        breaks = metadata[, SampleID],
        labels = metadata[, get("Patient ID")]
    )
    + scale_fill_manual(
        limits = metadata[, SampleID],
        labels = metadata[, get("Patient ID")],
        values = metadata[, Colour],
        name = "Patient"
    )
    + guides(fill = FALSE)
    + theme_minimal()
    + theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    )
)
savefig(gg_bounds_ctcf_fc, file.path(PLOT_DIR, "boundary-counts.ctcf-proximity.fold"))


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
    # + geom_path(aes(x = w, y = N, colour = Patient_ID, group = Patient_ID))
    + geom_point(
        aes(x = w, y = N, colour = Patient_ID),
        position = position_jitter(width = 0.2, height = 0)
    )
    + geom_boxplot(aes(x = w, y = N, group = w), outlier.shape = NA, alpha = 0.2)
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
    ggplot(data = tad_size_ecdf_est)
    + geom_path(aes(x = x, y = 100 * Mean, colour = w, group = w))
    + labs(x = "TAD Size (Mbp)", y = "% of TADs (Cumulative Density)")
    + scale_x_continuous(
        limits = c(0, 5e6),
        breaks = seq(0, 5e6, by = 1e6),
        labels = seq(0, 5)
    )
    + scale_colour_viridis_c(
        breaks = seq(MIN_WINDOW, MAX_WINDOW, length.out = 4),
        name = "Window Size"
    )
    + theme_minimal()
    + theme(
        legend.position = "bottom"
    )
)
savefig(gg_tad_size, file.path(PLOT_DIR, "tad-size.distribution"), height = 20)

# same as above but only focussing on a few window sizes
w_subset <- seq(MIN_WINDOW, MAX_WINDOW, length.out = 4)
gg_tad_size_reduced <- (
    ggplot(data = tad_size_ecdf_est[w %in% w_subset])
    + geom_path(
        aes(
            x = x,
            y = 100 * Mean,
            colour = w
        )
    )
    + labs(x = "TAD Size (Mbp)", y = "% of TADs (Cumulative Density)")
    + scale_x_continuous(
        limits = c(0, 5e6),
        breaks = seq(0, 5e6, by = 1e6),
        labels = seq(0, 5)
    )
    + scale_colour_viridis_c(
        breaks = seq(MIN_WINDOW, MAX_WINDOW, length.out = 4),
        name = "Window Size"
    )
    # + facet_wrap(~ w, nrow = 1)
    + theme_minimal()
    + theme(
        # legend.position = "bottom"
    )
)
savefig(gg_tad_size_reduced, file.path(PLOT_DIR, "tad-size.distribution.reduced.faceted"), height = 12, width = 30)

gg_tad_size_reduced <- (
    ggplot(data = tad_size_ecdf_est[w %in% w_subset])
    + geom_path(
        aes(
            x = x,
            y = 100 * Mean,
            colour = w,
            group = w
        )
    )
    + labs(x = "TAD Size (Mbp)", y = "% of TADs (Cumulative Density)")
    + scale_x_continuous(
        limits = c(0, 5e6),
        breaks = seq(0, 5e6, by = 1e6),
        labels = seq(0, 5)
    )
    + scale_colour_viridis_c(
        breaks = seq(MIN_WINDOW, MAX_WINDOW, length.out = 4),
        name = "Window Size"
    )
    + theme_minimal()
    + theme(
        # legend.position = "bottom"
    )
)
savefig(gg_tad_size_reduced, file.path(PLOT_DIR, "tad-size.distribution.reduced"), height = 14)

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
    + geom_boxplot(
        aes(x = w, y = 1 - dist, group = w),
        alpha = 0.5,
        outlier.shape = NA,
        width = 0.5
    )
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
