# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))

MAX_WINDOW = 28

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

# ==============================================================================
# Data
# ==============================================================================
# load metadata
metadata = fread(file.path("..", "..", "Data", "External", "LowC_Samples_Data_Available.tsv"))

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
tads[, lower_persistence := pmin(MAX_WINDOW, lower_persistence)]
tads[, upper_persistence := pmin(MAX_WINDOW, upper_persistence)]
boundaries[, Order := pmin(MAX_WINDOW, Order)]

# load BPscore calculations
bpscore <- fread("Statistics/tad-distances.tsv", sep = "\t", header = TRUE)
window_diffs <- fread("Statistics/tad-similarity-deltas.tsv", sep = "\t", header = TRUE)

# ==============================================================================
# Analysis
# ==============================================================================
boundary_counts = boundaries[, .N, by = "SampleID"]
boundary_counts_order = boundaries[, .N, by = c("SampleID", "Order")]

# calculate coefficient of variation across TAD sizes to see where samples vary
tads_mean_size <- tads[, mean(width), by = c("SampleID", "w")]
tads_cov <- tads_mean_size[, sd(V1) / mean(V1), by = "w"]
tads_median_size <- tads[, median(width), by = c("SampleID", "w")]
tads_madm <- tads_median_size[, median(abs(V1 - median(V1))), by = "w"]

# ==============================================================================
# Plots
# ==============================================================================
# plot number of resolved boundaries
gg_boundaries = (
    ggplot(data = boundary_counts)
    + geom_col(aes(x = SampleID, y = N, fill = SampleID))
    + scale_fill_viridis_d()
    + labs(x = NULL, y = "Number of unique boundaries")
    + guides(fill = FALSE)
    + theme_minimal()
    + theme(
        axis.text.x = element_text(angle = 90)
    )
)
ggsave(
    "Plots/tad-counts.png",
    gg_boundaries,
    height = 12,
    width = 20,
    units = "cm"
)
ggsave(
    "Plots/tad-counts.pdf",
    gg_boundaries,
    height = 12,
    width = 20,
    units = "cm",
    dpi = 400
)

# plot number of resolved boundaries by order
gg_bounds_order = (
    ggplot(data = boundary_counts_order)
    + geom_col(aes(x = Order, y = N, fill = SampleID, group = SampleID), position = "dodge")
    + scale_fill_viridis_d()
    + xlim(c(1, MAX_WINDOW))
    + labs(x = "Boundary Order", y = "Number of unique boundaries")
    + theme_minimal()
    + theme(
        axis.text.x = element_text(angle = 90)
    )
)
ggsave(
    "Plots/tad-counts-by-order.png",
    gg_bounds_order,
    height = 12,
    width = 20,
    units = "cm"
)
ggsave(
    "Plots/tad-counts-by-order.pdf",
    gg_bounds_order,
    height = 12,
    width = 20,
    units = "cm",
    dpi = 400
)


# unique TAD boundaries by the number of unique samples they appear in
gg_uniq_boundaries <- (
    ggplot(data = boundary_counts_recurrence[, .N, by = "Recurrence"])
    + geom_col(aes(x = Recurrence, y = N, fill = Recurrence))
    + geom_text(
        aes(
            x = Recurrence,
            y = N,
            label = paste0(N, "\n", "(", 100 * round(N / sum(N), 3), "%)")
        ),
        size = 2,
        vjust = -0.2
    )
    + labs(x = "Patients with the same boundary", y = "Frequency")
    + scale_fill_viridis_c()
    + scale_x_continuous(
        breaks = seq(1, 13, 3),
        labels = seq(1, 13, 3)
    )
    + scale_y_continuous(
        limits = c(0, 7100)
    )
    + guides(fill = FALSE)
    + theme_minimal()
)
ggsave(
    "Plots/tad-boundary-by-recurrence.png",
    gg_uniq_boundaries,
    height = 12,
    width = 20,
    units = "cm"
)


gg_tad_size = (
    ggplot(data = tads)
    #+ geom_density(aes(x = (end - start) / 1e6, colour = SampleID), alpha = 0.1)
    + geom_density(aes(x = width / 1e6, colour = SampleID), alpha = 0.1)
    #+ geom_smooth(aes(x = width / 1e6, y = N))
    + labs(x = "TAD Size (Mbp)", y = "Scaled density")
    + xlim(c(0, 5))
    + facet_wrap(~ w, scales = "free_y", nrow = 5)
    + theme_minimal()
    + theme(
        legend.position = "bottom"
    )
)
ggsave(
    "Plots/tad-size-distribution.png",
    gg_tad_size,
    height = 20,
    width = 20,
    units = "cm"
)
ggsave(
    "Plots/tad-size-distribution.pdf",
    gg_tad_size,
    height = 20,
    width = 20,
    units = "cm",
    dpi = 400
)

gg = (
    ggplot(data = tads_cov)
    + geom_col(aes(x = w, y = V1))
    + labs(x = "Window size", y = "Coefficient of Variation")
    + theme_minimal()
)
ggsave(
    "Plots/tad-size.coeff-var.png",
    height = 20,
    width = 20,
    units = "cm"
)

# plot distances between TADs of different samples across window sizes
gg_bpscore = (
    ggplot(data = bpscore[w < MAX_WINDOW])
    + geom_point(aes(x = w, y = 1 - dist, colour = w), position = position_jitter(width = 0.2, height = 0))
    #+ geom_violin(aes(x = w, y = dist, group = w, fill = w))
    + geom_boxplot(aes(x = w, y = 1 - dist, group = w), alpha = 0.5, outlier.shape = NA, width = 0.5)
    + ylim(c(0.68, 1))
    + labs(x = "Window size", y = "TAD similarity (1 - BPscore)")
    + scale_colour_viridis_c()
    + scale_fill_viridis_c()
    + guides(colour = FALSE)
    + theme_minimal()
)
ggsave(
    "Plots/bp-score.png",
    gg_bpscore,
    height = 12,
    width = 20,
    units = "cm"
)
ggsave(
    "Plots/bp-score.pdf",
    gg_bpscore,
    height = 12,
    width = 20,
    units = "cm",
    dpi = 400
)


# plot the change in similarities over window sizes
gg_sim_window = (
    ggplot(data = window_diffs)
    + geom_path(aes(x = w, y = 1 - diff, group = SampleID, colour = SampleID))
    + labs(x = "Window size", y = expression("1 - " * delta[w]))
    + theme_minimal()
)
ggsave(
    "Plots/bp-score.window-similarity.png",
    gg_sim_window,
    height = 12,
    width = 20,
    units = "cm"
)
ggsave(
    "Plots/bp-score.window-similarity.pdf",
    gg_sim_window,
    height = 12,
    width = 20,
    units = "cm",
    dpi = 400
)


# plot the change in similarities over window sizes
gg_sim_window_delta = (
    ggplot(data = window_diffs[w >= 5])
    + geom_path(aes(x = w, y = abs_delta, group = SampleID, colour = SampleID))
    + labs(x = "Window size", y = expression("|" * delta[w] - delta[w-1] * "|"))
    + theme_minimal()
)
ggsave(
    "Plots/bp-score.window-similarity.delta.png",
    gg_sim_window_delta,
    height = 12,
    width = 20,
    units = "cm"
)
ggsave(
    "Plots/bp-score.window-similarity.delta.pdf",
    gg_sim_window_delta,
    height = 12,
    width = 20,
    units = "cm",
    dpi = 400
)
