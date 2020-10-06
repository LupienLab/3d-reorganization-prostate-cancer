# ==============================================================================
# Meta
# ==============================================================================
# Get loop anchors
# --------------------------------------
# Description: Take 2D loop calls and identify the unique anchors
# Author: James Hawley


# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("logging"))

loginfo("Loading packages")
suppressMessages(library("data.table"))
suppressMessages(library("GenomicRanges"))
suppressMessages(library("regioneR"))
suppressMessages(library("ggplot2"))
suppressMessages(library("ggupset"))

# ==============================================================================
# Functions
# ==============================================================================
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
loginfo("Loading data")
# load metadata
metadata <- fread("config.tsv", sep = "\t", header = TRUE)
metadata <- metadata[Include == "Yes", .SD]
SAMPLES <- metadata[, SampleID]

# load loop calls
loops <- rbindlist(lapply(
    SAMPLES,
    function(s) {
        dt <- fread(
            file.path("..", "..", "Data", "Processed", "2019-06-18_PCa-LowC-sequencing", "Loops", paste0(s, ".res_10000bp.loops.tsv")),
            sep = "\t",
            header = TRUE
        )
        dt[, SampleID := s]
        return(dt)
    }
))

loops_gr <- lapply(
    SAMPLES,
    function(s) {
        # get loops just for a single sample
        dt <- loops[SampleID == s, .SD]
        # convert to GRanges object
    }
)

# convert to GRanges for overlap calculations
loop_anchors <- reduce(
    GRanges(
        seqnames = loops[, c(chr_x, chr_y)],
        ranges = IRanges(
            start = loops[, c(start_x, start_y)],
            end = loops[, c(end_x, end_y)]
        )
    )
)

# ==============================================================================
# Analysis
# ==============================================================================
loginfo("Merging loop calls")

# map original loop anchors to reduced anchors
loops[, `:=`(
    # get subjectHits (from `loop_anchors`)
    anchor_ID_x = subjectHits(findOverlaps(
        # create GRanges object for overlaps
        GRanges(
            seqnames = chr_x,
            ranges = IRanges(
                start = start_x,
                end = end_x
            )
        ),
        # overlap with loop anchors
        loop_anchors
    )),
    # same as above but with the second anchor in the loop
    anchor_ID_y = subjectHits(findOverlaps(
        GRanges(
            seqnames = chr_y,
            ranges = IRanges(
                start = start_y,
                end = end_y
            )
        ),
        loop_anchors
    ))
)]

# concatenate columns to create a unique loop ID
loops[, loopID := paste0(anchor_ID_x, "_", anchor_ID_y)]

# merge resulting loop calls to give a set of loops based on the anchors from all samples
corrected_anchors_x <- loop_anchors[loops[, anchor_ID_x]]
corrected_anchors_y <- loop_anchors[loops[, anchor_ID_y]]

merged_loops <- data.table(
    "SampleID" = loops[, SampleID],
    "chr_x" = as.vector(seqnames(corrected_anchors_x)),
    "start_x" = start(corrected_anchors_x),
    "end_x" = end(corrected_anchors_x),
    "chr_y" = as.vector(seqnames(corrected_anchors_y)),
    "start_y" = start(corrected_anchors_y),
    "end_y" = end(corrected_anchors_y),
    "anchor_ID_x" = loops[, anchor_ID_x],
    "anchor_ID_y" = loops[, anchor_ID_y],
    "loopID" = loops[, loopID],
    "detection_scale" = loops[, detection_scale],
    "fdr" = loops[, fdr]
)

# ==============================================================================
# Save data
# ==============================================================================
loginfo("Saving data")
fwrite(merged_loops, file.path("Loops", "merged-loops.tsv"), sep = "\t", col.names = TRUE)

# ==============================================================================
# Plots
# ==============================================================================
loginfo("Plotting figures")

# collapse list of loopIDs by the sample
collapsed_IDs_by_sample <- merged_loops[, .(SampleID = list(unique(SampleID))), keyby = "loopID"]

# collapse list of loopIDs by the tissue type
collapsed_IDs_by_type <- merge(
    merged_loops,
    metadata[, .SD, .SDcols = c("SampleID", "Type")]
)[, .(Type = list(unique(Type))), keyby = "loopID"]

# plots
gg_sample <- (
    ggplot(data = collapsed_IDs_by_sample)
    + geom_bar(aes(x = SampleID))
    + labs(x = "Samples", y = "Shared Loops")
    + scale_x_upset(order_by = "freq", n_intersections = Inf)
    + theme_minimal()
)
savefig(gg_sample, file.path("Plots", "loop-calls.by-sample"))

gg_type <- (
    ggplot(data = collapsed_IDs_by_type)
    + geom_bar(aes(x = Type))
    + labs(x = "Tissue", y = "Shared Loops")
    + scale_x_upset(order_by = "freq", n_intersections = Inf)
    + theme_minimal()
)
savefig(gg_type, file.path("Plots", "loop-calls.by-type"))
