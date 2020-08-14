# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
suppressMessages(library("scales"))
suppressMessages(library("MASS"))
suppressMessages(library("pheatmap"))

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
cat("Loading Data\n")
# load metadata
metadata <- fread("config.tsv", sep = "\t", header = TRUE)
metadata <- metadata[Include == "Yes"]
SAMPLES <- metadata[, SampleID]
TUMOUR_SAMPLES <- metadata[Source == "Primary" & Type == "Malignant", SampleID]
LINE_SAMPLES <- metadata[Source == "Cell Line", SampleID]

# load aggregated boundary calls from each sample
boundaries = rbindlist(lapply(
    TUMOUR_SAMPLES,
    function(s) {
        dt = fread(
            file.path("resolved-TADs", paste0(s, ".40000bp.aggregated-boundaries.tsv")),
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

# only keep TADs called at w <= MAX_WINDOW and no chrY
tads <- tads[w <= MAX_WINDOW & chr != "chrY"]
tads[, lower_persistence := pmin(MAX_PERSISTENCE, lower_persistence)]
tads[, upper_persistence := pmin(MAX_PERSISTENCE, upper_persistence)]
boundaries[, Persistence := pmin(MAX_PERSISTENCE, Order)]

# add Patient IDs to tads and boundaries, not just SampleIDs
tads <- merge(tads, metadata[, .SD, .SDcols = c("SampleID", "Label", "Type_Colour", "Sample_Colour")])
boundaries <- merge(boundaries, metadata[, .SD, .SDcols = c("SampleID", "Label", "Type_Colour", "Sample_Colour")])

# load BPscore calculations
bpscore <- fread("Statistics/tad-distances.tsv", sep = "\t", header = TRUE)
window_diffs <- fread("Statistics/tad-similarity-deltas.tsv", sep = "\t", header = TRUE)

# TAD summary
tad_summary <- tads[, .N, keyby = c("SampleID", "w")]
tad_summary <- merge(tad_summary, metadata[, .SD, .SDcols = c("SampleID", "Label")], by = "SampleID")
fwrite(tad_summary, "Statistics/tad-counts.tsv", sep = "\t", col.names = TRUE)

# ==============================================================================
# Analysis
# ==============================================================================
cat("Counting boundaries\n")
boundary_counts = boundaries[, .N, by = c("SampleID", "Label")]
boundary_counts_persistence = boundaries[, .N, by = c("SampleID", "Label", "Persistence")]
boundary_max_persistence <- merge(
    boundaries[Persistence == MAX_PERSISTENCE, .(Max_Persistence = .N), by = "SampleID"],
    boundaries[Persistence < MAX_PERSISTENCE, .(Lesser_Persistence = .N), by = "SampleID"]
)
boundary_max_persistence[, Frac_Max_Persistence := Max_Persistence / (Max_Persistence + Lesser_Persistence)]
fwrite(boundary_max_persistence, "Statistics/boundary-hierarchy.tsv", sep = "\t", col.names = TRUE)

cat("Counting TADs\n")
# calculate coefficient of variation across TAD sizes to see where samples vary
tads_mean_size <- tads[, mean(width), by = c("SampleID", "w")]
tads_cov <- tads_mean_size[, sd(V1) / mean(V1), by = "w"]
tads_median_size <- tads[, median(width), by = c("SampleID", "w")]
tads_madm <- tads_median_size[, median(abs(V1 - median(V1))), by = "w"]

# ECDFs number of TADs of a given width
tad_size_ecdf <- rbindlist(lapply(
    TUMOUR_SAMPLES,
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

# colour bpscore rows by primary-primary, line-line, or primary-line comparisons
cat("Counting similarity\n")
bpscore <- merge(
    bpscore,
    metadata[, .(SampleID, Source, Label, Type_Colour)],
    by.x = "s2",
    by.y = "SampleID",
    all.x = FALSE,
    all.y = TRUE
)
bpscore <- merge(
    bpscore,
    metadata[, .(SampleID, Source, Label, Type_Colour)],
    by.x = "s1",
    by.y = "SampleID",
    suffixes = c("_2", "_1"),
    all.x = FALSE,
    all.y = TRUE,
)
# if s1 and s2 come from the same source material, return Type_Colour, otherwise return green colour
bpscore[Source_1 == Source_2 & Source_1 == "Primary", Colour := "#1F77B4"]
bpscore[Source_1 == Source_2 & Source_1 == "Cell Line", Colour := "#FF7F0D"]
bpscore[Source_1 != Source_2, Colour := "#3FE686"]

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
       labels = metadata[, Label]
   )
   + scale_fill_manual(
       limits = metadata[, SampleID],
       values = metadata[, Type_Colour]
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
        labels = metadata[, Label],
        values = metadata[, Type_Colour],
        name = "Patient"
    )
    + theme_minimal()
    + theme(
        axis.text.x = element_text(angle = 90)
    )
)
savefig(gg_bounds_persistence, file.path(PLOT_DIR, "boundary-counts.by-persistence"))

# TADs
# --------------------------------------
# TAD counts per sample faceted by window size
gg_tad_counts_window <- (
    ggplot(data = tads[, .N, by = c("SampleID", "w", "Sample_Colour")])
    + geom_col(aes(x = SampleID, y = N, fill = Sample_Colour))
    + labs(x = NULL, y = "Number of TADs")
    + scale_fill_manual(
        limits = metadata[, Sample_Colour],
        labels = metadata[, paste(Source, Type)],
        values = metadata[, Sample_Colour]
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
    ggplot(data = tads[, .N, by = c("SampleID", "w")])
    + geom_col(aes(x = w, y = N, fill = w))
    + labs(x = "Window size", y = "Number of TADs")
    + scale_fill_viridis_c()
    + guides(fill = FALSE)
    + facet_wrap(
        ~ SampleID,
        nrow = 3
        # labeller = as_labeller(metadata[, Label])
    )
    + theme_minimal()
    + theme(
        legend.position = "bottom"
    )
)
savefig(gg_tad_counts_sample, file.path(PLOT_DIR, "tad-counts.by-sample"), height = 20)

# TAD counts per sample at each window size
gg_tad_counts <- (
    ggplot(data = tads[SampleID %in% TUMOUR_SAMPLES, .N, by = c("SampleID", "w", "Sample_Colour")])
    + geom_point(
        aes(x = w, y = N, colour = SampleID),
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
        limits = metadata[SampleID %in% TUMOUR_SAMPLES, SampleID],
        labels = metadata[SampleID %in% TUMOUR_SAMPLES, Label],
        values = metadata[SampleID %in% TUMOUR_SAMPLES, Sample_Colour],
        name = "Source"
    )
    + theme_minimal()
)
savefig(gg_tad_counts, file.path(PLOT_DIR, "tad-counts"))

# TAD size distribution
ribbon_alpha <- 0.4
gg_tad_size <- (
    ggplot(data = tad_size_ecdf_est)
    + geom_path(aes(x = x, y = 100 * Mean, colour = w, group = w))
    + geom_ribbon(
        aes(x = x, ymin = 100 * (Mean - SD), ymax = 100 * (Mean + SD), fill = w, group = w),
        alpha = ribbon_alpha
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
    + scale_fill_viridis_c(
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
    + geom_path(aes(x = x, y = 100 * Mean, colour = w, group = w))
    + geom_ribbon(
        aes(x = x, ymin = 100 * (Mean - SD), ymax = 100 * (Mean + SD), fill = w, group = w),
        alpha = ribbon_alpha
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
    + scale_fill_viridis_c(
        breaks = seq(MIN_WINDOW, MAX_WINDOW, length.out = 4),
        name = "Window Size"
    )
    + theme_minimal()
    + theme(
        legend.position = "bottom"
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
    ggplot(data = bpscore[s1 != s2 & w <= MAX_WINDOW & !grepl("(4DN|SRR|BP)", s1) & !grepl("(4DN|SRR|BP)", s2)])
    + geom_point(
        # aes(x = w, y = 1 - dist, colour = Colour),
        aes(x = w, y = 1 - dist, colour = w),
        position = position_jitter(width = 0.2, height = 0)
    )
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
    # + scale_colour_manual(
    #     limits = bpscore[, unique(Colour)],
    #     values = bpscore[, unique(Colour)]
    # )
    + scale_colour_viridis_c()
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
        labels = metadata[, Label],
        values = metadata[, Sample_Colour],
        name = "Sample"
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
        labels = metadata[, Label],
        values = metadata[, Sample_Colour],
        name = "Patient"
    )
    + theme_minimal()
)
savefig(gg_sim_window_delta, file.path(PLOT_DIR, "bp-score.window-similarity.delta"))

# cluster samples by BPscore distances
# recast data to wide form for plotting with pheatmap
mat <- as.matrix(
    dcast(bpscore[w == MAX_WINDOW], s1 ~ s2, value.var = "dist", fill = 0),
    rownames = "s1"
)
annot_col <- as.data.frame(metadata[, .SD, .SDcols = c("Type", "Tissue", "Source")])
rownames(annot_col) <- metadata[, SampleID]

for (ext in c("png", "pdf")) {
    pheatmap(
        mat = 1 - mat,
        filename = paste0(file.path(PLOT_DIR, "bp-score.cluster."), ext),
        clustering_rows = TRUE,
        clustering_cols = TRUE,
        clustering_distance_rows = as.dist(mat),
        clustering_distance_cols = as.dist(mat),
        clustering_method = "ward.D2",
        labels_row = metadata[order(SampleID), Label],
        labels_col = metadata[order(SampleID), Label],
        legend = TRUE,
        show_rownames = FALSE,
        annotation_col = annot_col
    )
}
