# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("scales"))
suppressMessages(library("ggplot2"))

MAX_WINDOW <- 20
MIN_WINDOW <- 3
RESOLUTION <- 40000
MAX_PERSISTENCE <- MAX_WINDOW - MIN_WINDOW + 1
PLOT_DIR <- "Plots"

# ==============================================================================
# Functions
# ==============================================================================
split_comma_col <- function(v, f=identity) {
    # split into list
    splitv <- lapply(v, function(x) {strsplit(x, ",")[[1]]})
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
savefig <- function(gg, prefix, ext = c("png", "pdf"), width = 20, height = 12, dpi = 400) {
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

# load aggregated boundary calls from each sample
boundaries <- rbindlist(lapply(
    PRIMARY_SAMPLES,
    function(s) {
        dt = fread(
            file.path("Aggregated-TADs", paste0(s, ".300000000.res_40000bp.agg-boundaries.tsv")),
            sep = "\t",
            header = TRUE
        )
        dt[, SampleID := s]
        return(dt)
    }
))
boundaries$w <- split_comma_col(boundaries$w, as.numeric)
boundaries$left_type <- split_comma_col(boundaries$left_type)
boundaries$right_type <- split_comma_col(boundaries$right_type)

# merge metadata information into boundaries
boundaries <- merge(
    x = boundaries,
    y = metadata[, .SD, .SDcols = c("SampleID", "Label", "Type")],
    by = "SampleID"
)

# ==============================================================================
# Analysis
# ==============================================================================
cat("Counting boundaries\n")
boundary_counts <- boundaries[, .N, by = c("SampleID", "Label", "Type")]
boundary_counts_persistence = boundaries[, .N, by = c("SampleID", "Label", "Type", "Persistence")]
boundary_max_persistence <- merge(
    boundaries[Persistence == MAX_PERSISTENCE, .(Max_Persistence = .N), by = "SampleID"],
    boundaries[Persistence < MAX_PERSISTENCE, .(Lesser_Persistence = .N), by = "SampleID"]
)
boundary_max_persistence[, Frac_Max_Persistence := Max_Persistence / (Max_Persistence + Lesser_Persistence)]
fwrite(boundary_max_persistence, "Statistics/boundary-hierarchy.tsv", sep = "\t", col.names = TRUE)

boundaries_singleton <- boundaries[Persistence == 1]

# ==============================================================================
# Plots
# ==============================================================================
# Boundaries
# --------------------------------------
# plot number of resolved boundaries
gg_boundaries <- (
    ggplot(data = boundary_counts)
    + geom_col(aes(x = SampleID, y = N, fill = Type))
    + scale_x_discrete(
        limits = metadata[, SampleID],
        labels = metadata[, Label]
    )
    + scale_fill_manual(
        limits = c("Benign", "Malignant"),
        labels = c("Benign", "Tumour"),
        values = c("#AEC7E8", "#1F77B4"),
        name = "Sample Type"
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
    + geom_col(aes(x = Persistence, y = N, fill = Type, group = SampleID), position = "dodge")
    + labs(x = "Boundary Persistence", y = "Number of unique boundaries")
    + scale_x_discrete(
        breaks = seq(1, MAX_PERSISTENCE, length.out = 4),
        limits = seq(1, MAX_PERSISTENCE, length.out = 4),
        labels = seq(1, MAX_PERSISTENCE, length.out = 4)
    )
    + scale_fill_manual(
        limits = c("Benign", "Malignant"),
        labels = c("Benign", "Tumour"),
        values = c("#AEC7E8", "#1F77B4"),
        name = "Sample Type"
    )
    + theme_minimal()
    + theme(
        axis.text.x = element_text(angle = 90),
        legend.position = "bottom"
    )
)
savefig(gg_bounds_persistence, file.path(PLOT_DIR, "boundary-counts.by-persistence"))
