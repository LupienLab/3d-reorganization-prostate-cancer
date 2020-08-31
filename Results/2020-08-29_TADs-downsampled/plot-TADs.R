# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("scales"))
suppressMessages(library("ggplot2"))

MAX_WINDOW <- 23
MIN_WINDOW <- 3
RESOLUTION <- 40000
MAX_PERSISTENCE <- MAX_WINDOW - MIN_WINDOW + 1
PLOT_DIR <- "Plots"

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
cat("Loading Data\n")
# load metadata
metadata <- fread("config.tsv", sep = "\t", header = TRUE)
metadata <- metadata[Include == "Yes"]
SAMPLES <- metadata[, SampleID]
TUMOUR_SAMPLES <- metadata[Source == "Primary" & Type == "Malignant", SampleID]
BENIGN_SAMPLES <- metadata[Source == "Primary" & Type == "Benign", SampleID]
PRIMARY_SAMPLES <- metadata[Source == "Primary", SampleID]

# load TADs
tads <- rbindlist(lapply(
    SAMPLES,
    function(s) {
        dt <- fread(
            file.path("Aggregated-TADs", paste0(s, ".300000000.res_", RESOLUTION, "bp.agg-domains.tsv")),
            sep = "\t",
            header = TRUE
        )
        dt[, SampleID := s]
        return(dt)
    }
))
tads[, width := as.numeric(end - start)]

# add labels to TADs and boundaries, not just SampleIDs
tads <- merge(tads, metadata[, .SD, .SDcols = c("SampleID", "Label", "Type_Colour", "Sample_Colour", "Type")])
# boundaries <- merge(boundaries, metadata[, .SD, .SDcols = c("SampleID", "Label", "Type_Colour", "Sample_Colour")])

# TAD summary (grepl => only keep TADs labelled domains, no gaps)
tad_summary <- tads[grepl("domain", type), .(N_TADs = .N), keyby = c("w", "SampleID", "Label", "Type")]
fwrite(tad_summary, "Statistics/tad-counts.tsv", sep = "\t", col.names = TRUE)

# ==============================================================================
# Analysis
# ==============================================================================
cat("Counting TADs\n")
# calculate coefficient of variation across TAD sizes to see where samples vary
tads_mean_size <- tads[, mean(width), by = c("SampleID", "w")]
tads_cov <- tads_mean_size[, sd(V1) / mean(V1), by = "w"]
tads_median_size <- tads[, median(width), by = c("SampleID", "w")]
tads_madm <- tads_median_size[, median(abs(V1 - median(V1))), by = "w"]

# ECDFs number of TADs of a given width
tad_size_ecdf <- rbindlist(lapply(
    PRIMARY_SAMPLES,
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


# ==============================================================================
# Plots
# ==============================================================================
# TAD counts per sample faceted by window size
gg_tad_counts_window <- (
    ggplot(data = tad_summary)
    + geom_col(aes(x = SampleID, y = N_TADs, fill = SampleID))
    + labs(x = NULL, y = "Number of TADs")
    + scale_x_discrete(
        limits = metadata[, SampleID],
        labels = metadata[, Label]
    )
    + scale_fill_manual(
        limits = metadata[, SampleID],
        labels = metadata[, paste(Source, Type)],
        values = metadata[, Type_Colour]
    )
    + guides(fill = FALSE)
    + facet_wrap(~ w, nrow = 6)
    + theme_minimal()
    + theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "bottom"
    )
)
savefig(gg_tad_counts_window, file.path(PLOT_DIR, "tad-counts.by-window"), height = 30, width = 30)

# TAD counts per sample faceted by window size
gg_tad_counts_sample <- (
    ggplot(data = tad_summary)
    + geom_col(aes(x = w, y = N_TADs, fill = w))
    + labs(x = "Window size", y = "Number of TADs")
    + scale_fill_viridis_c()
    + guides(fill = FALSE)
    + facet_wrap(
        ~ Label,
        nrow = 3
        # labeller = as_labeller(metadata[, Label])
    )
    + theme_minimal()
    + theme(
        legend.position = "bottom"
    )
)
savefig(gg_tad_counts_sample, file.path(PLOT_DIR, "tad-counts.by-sample"), height = 20, width = 30)

# TAD counts per sample at each window size
gg_tad_counts <- (
    ggplot(
        data = tad_summary,
        mapping = aes(
            x = w,
            y = N_TADs,
            colour = Type,
            group = paste(w,  Type)
        )
    )
    + geom_point(
        position = position_jitterdodge(
            jitter.width = 0.2, jitter.height = 0,
            dodge.width = 0.9
        ),
        size = 2
    )
    + geom_boxplot(
        position = position_dodge(),
        outlier.shape = NA,
        alpha = 0.3
    )
    + labs(x = "Window size", y = "Number of TADs")
    + scale_x_discrete(
        breaks = seq(MIN_WINDOW, MAX_WINDOW, length.out = 5),
        limits = seq(MIN_WINDOW, MAX_WINDOW, length.out = 5),
        labels = seq(MIN_WINDOW, MAX_WINDOW, length.out = 5)
    )
    + scale_colour_manual(
        limits = c("Benign", "Malignant"),
        labels = c("Benign", "Tumour"),
        values = c("#AEC7E8", "#1F77B4"),
        name = "Sample Type"
    )
    + theme_minimal()
    + theme(legend.position = "bottom")
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

# same as above but only focusing on a few window sizes
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
