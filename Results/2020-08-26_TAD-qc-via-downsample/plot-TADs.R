# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("scales"))
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
cat("Loading Data\n")
# load metadata
metadata <- fread("config.tsv", sep = "\t", header = TRUE)
metadata <- metadata[Include == "Yes"]
SAMPLES <- metadata[, SampleID]

# load aggregated boundary calls from each sample
boundaries = rbindlist(lapply(
    SAMPLES,
    function(s) {
        dt = fread(
            file.path("TADs", "resolved-TADs", paste0(s, ".40000bp.aggregated-boundaries.tsv")),
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
        dt <- fread(paste0("TADs/resolved-TADs/", s, ".40000bp.aggregated-domains.bedGraph"),
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

# TAD summary
tad_summary <- tads[, .N, keyby = c("SampleID", "w")]
tad_summary <- merge(tad_summary, metadata[, .SD, .SDcols = c("SampleID", "Label")], by = "SampleID")
fwrite(tad_summary, "Statistics/tad-counts.tsv", sep = "\t", col.names = TRUE)

# ==============================================================================
# Analysis
# ==============================================================================
cat("Counting boundaries\n")
boundary_counts = boundaries[, .N, by = c("SampleID", "Label", "Type_Colour")]
boundary_counts_persistence = boundaries[, .N, by = c("SampleID", "Label", "Persistence", "Type_Colour")]
boundary_max_persistence <- merge(
    boundaries[Persistence == MAX_PERSISTENCE, .(Max_Persistence = .N), by = "SampleID"],
    boundaries[Persistence < MAX_PERSISTENCE, .(Lesser_Persistence = .N), by = "SampleID"]
)
boundary_max_persistence[, Frac_Max_Persistence := Max_Persistence / (Max_Persistence + Lesser_Persistence)]
fwrite(boundary_max_persistence, "Statistics/boundary-hierarchy.tsv", sep = "\t", col.names = TRUE)

boundaries_singleton <- boundaries[Persistence == 1]

cat("Counting TADs\n")
# calculate coefficient of variation across TAD sizes to see where samples vary
tads_mean_size <- tads[, mean(width), by = c("SampleID", "w")]
tads_cov <- tads_mean_size[, sd(V1) / mean(V1), by = "w"]
tads_median_size <- tads[, median(width), by = c("SampleID", "w")]
tads_madm <- tads_median_size[, median(abs(V1 - median(V1))), by = "w"]


# ==============================================================================
# Plots
# ==============================================================================
# Boundaries
# --------------------------------------
# plot number of resolved boundaries
gg_boundaries <- (
    ggplot(data = boundary_counts)
    + geom_col(aes(x = SampleID, y = N, fill = Type_Colour))
   + scale_x_discrete(
       limits = metadata[, SampleID],
       labels = metadata[, Label]
   )
   + scale_fill_manual(
       limits = c("#1F77B4", "#AEC7E8"),
       labels = c("Tumour", "Benign"),
       values = c("#1F77B4", "#AEC7E8"),
       name = "Source"
   )
    + labs(x = NULL, y = "Number of unique boundaries")
    + theme_minimal()
    + theme(
        axis.text.x = element_text(angle = 90),
        legend.position = "bottom"
    )
)
savefig(gg_boundaries, file.path(PLOT_DIR, "boundary-counts"))

# plot number of resolved boundaries by order
gg_bounds_persistence <- (
    ggplot(data = boundary_counts_persistence)
    + geom_col(aes(x = Persistence, y = N, fill = Type_Colour, group = SampleID), position = "dodge")
    + labs(x = "Boundary Persistence", y = "Number of unique boundaries")
    + scale_x_discrete(
        breaks = c(1, 6, 11, 16, MAX_PERSISTENCE),
        limits = c(1, 6, 11, 16, MAX_PERSISTENCE),
        labels = c(1, 6, 11, 16, MAX_PERSISTENCE)
    )
    + scale_fill_manual(
        limits = c("#1F77B4", "#AEC7E8"),
        labels = c("Tumour", "Benign"),
        values = c("#1F77B4", "#AEC7E8"),
        name = "Source"
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
    ggplot(data = tads[, .N, by = c("SampleID", "w", "Type_Colour")])
    + geom_col(aes(x = SampleID, y = N, fill = Type_Colour))
    + labs(x = NULL, y = "Number of TADs")
    + scale_x_discrete(
        limits = metadata[, SampleID],
        labels = metadata[, Label]
    )
    + scale_fill_manual(
        limits = metadata[, Type_Colour],
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
    ggplot(data = tads[, .N, by = c("SampleID", "Label", "w")])
    + geom_col(aes(x = w, y = N, fill = w))
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
        data = tads[,
            .N,
            keyby = c("SampleID", "w", "Type_Colour")
        ],
        mapping = aes(
            x = w,
            y = N,
            colour = Type_Colour,
            group = paste(w, Type_Colour)
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
        breaks = seq(MIN_WINDOW, MAX_WINDOW, by = 3),
        limits = seq(MIN_WINDOW, MAX_WINDOW, by = 3),
        labels = seq(MIN_WINDOW, MAX_WINDOW, by = 3)
    )
    + scale_colour_manual(
        limits = c("#AEC7E8", "#1F77B4"),
        labels = c("Benign", "Tumour"),
        values = c("#AEC7E8", "#1F77B4"),
        name = "Sample Type"
    )
    + theme_minimal()
    + theme(legend.position = "bottom")
)
savefig(gg_tad_counts, file.path(PLOT_DIR, "tad-counts"))

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
