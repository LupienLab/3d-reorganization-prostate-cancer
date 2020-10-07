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


# ==============================================================================
# Analysis
# ==============================================================================
loginfo("Merging loop calls")

# calculate the resulting overlapping loop anchors
loop_anchors <- reduce(
    GRanges(
        seqnames = loops[, c(chr_x, chr_y)],
        ranges = IRanges(
            start = loops[, c(start_x, start_y)],
            end = loops[, c(end_x, end_y)]
        )
    )
)
# add loopID as a metadata column
mcols(loop_anchors)$anchorID <- 1:length(loop_anchors)

# map original loop anchors to reduced anchors
hits_x <- findOverlaps(
    # create GRanges object for overlaps
    GRanges(
        seqnames = loops[, chr_x],
        ranges = IRanges(
            start = loops[, start_x],
            end = loops[, end_x]
        )
    ),
    # overlap with loop anchors
    loop_anchors
)
hits_y <- findOverlaps(
    # create GRanges object for overlaps
    GRanges(
        seqnames = loops[, chr_y],
        ranges = IRanges(
            start = loops[, start_y],
            end = loops[, end_y]
        )
    ),
    # overlap with loop anchors
    loop_anchors
)

# assign the anchorIDs to the loop calls
loops[, `:=`(
    anchorID_x = subjectHits(hits_x),
    anchorID_y = subjectHits(hits_y),
    # concatenate columns to create a unique loop ID
    loopID = paste0(subjectHits(hits_x), "_", subjectHits(hits_y))
)]

# merge in tissue type information into the loops
loops <- merge(
    x = loops,
    y = metadata[, .SD, .SDcols = c("SampleID", "Type")],
    by = "SampleID"
)

# merge resulting loop calls to give a set of loops based on the anchors from all samples
corrected_anchors_x <- loop_anchors[loops[, anchorID_x]]
corrected_anchors_y <- loop_anchors[loops[, anchorID_y]]

merged_loops <- data.table(
    "SampleID" = loops[, SampleID],
    "Type" = loops[, Type],
    "chr_x" = as.vector(seqnames(corrected_anchors_x)),
    "start_x" = start(corrected_anchors_x),
    "end_x" = end(corrected_anchors_x),
    "chr_y" = as.vector(seqnames(corrected_anchors_y)),
    "start_y" = start(corrected_anchors_y),
    "end_y" = end(corrected_anchors_y),
    "anchor_ID_x" = loops[, anchorID_x],
    "anchor_ID_y" = loops[, anchorID_y],
    "loopID" = loops[, loopID],
    "detection_scale" = loops[, detection_scale],
    "fdr" = loops[, fdr]
)

# count the number of times a loop occurs
loop_counts <- dcast(
    # count loops byt their occurrences in each tissue type
    merged_loops[, .N, keyby = c("loopID", "Type")],
    # pivot the table to a wide format, rows are loops, columns are counts in benign and tumour
    loopID ~ Type,
    value.var = "N",
    # fill any NAs as 0
    fill = 0
)

# calculate if a loop occurs in at least 3 samples in one tissue type
loop_counts[, Consistent_in_benign_or_tumour := ((Benign >= 3) | (Malignant >= 3))]

# merge back in position information for loops
loop_counts_with_pos <- merge(
    x = loop_counts,
    y = merged_loops[, unique(.SD), .SDcols = c("chr_x", "start_x", "end_x", "chr_y", "start_y", "end_y", "loopID")],
    by = "loopID",
    all.x = TRUE,
    all.y = FALSE
)

# merge consistent calling in one tissue type back into merged loop calls
merged_loops_with_consistency <- merge(
    x = merged_loops,
    y = loop_counts,
    by = "loopID"
)


# collapse list of loopIDs by the sample
collapsed_IDs_by_sample <- merged_loops_with_consistency[
    Consistent_in_benign_or_tumour == TRUE,
    .(SampleID = list(unique(SampleID))),
    keyby = "loopID"
]

# collapse list of loopIDs by the tissue type
collapsed_IDs_by_type <- merged_loops_with_consistency[
    Consistent_in_benign_or_tumour == TRUE,
    .(Type = list(unique(Type))),
    keyby = "loopID"
]

# ==============================================================================
# Save data
# ==============================================================================
loginfo("Saving data")
fwrite(
    merged_loops_with_consistency,
    file.path("Loops", "merged-loops.tsv"),
    sep = "\t",
    col.names = TRUE
)
fwrite(
    # put the columns in an order compatible with BEDPE format
    loop_counts_with_pos[, .SD, .SDcols = c("chr_x", "start_x", "end_x", "chr_y", "start_y", "end_y", "loopID", "Benign", "Malignant", "Consistent_in_benign_or_tumour")],
    file.path("Loops", "merged-loops.sample-counts.tsv"),
    sep = "\t",
    col.names = TRUE
)

# save benign-specific loops
fwrite(
    # put the columns in an order compatible with BEDPE format
    loop_counts_with_pos[
        Consistent_in_benign_or_tumour == TRUE & Benign > 0 & Malignant == 0,
        .SD,
        .SDcols = c("chr_x", "start_x", "end_x", "chr_y", "start_y", "end_y", "loopID")
    ],
    file.path("Loops", "benign-specific-loops.tsv"),
    sep = "\t",
    col.names = TRUE
)

# save tumour-specific loops
fwrite(
    # put the columns in an order compatible with BEDPE format
    loop_counts_with_pos[
        Consistent_in_benign_or_tumour == TRUE & Benign == 0 & Malignant > 0,
        .SD,
        .SDcols = c("chr_x", "start_x", "end_x", "chr_y", "start_y", "end_y", "loopID")
    ],
    file.path("Loops", "tumour-specific-loops.tsv"),
    sep = "\t",
    col.names = TRUE
)

# save shared tumour-benign loops
fwrite(
    # put the columns in an order compatible with BEDPE format
    loop_counts_with_pos[
        Consistent_in_benign_or_tumour == TRUE & Benign > 0 & Malignant > 0,
        .SD,
        .SDcols = c("chr_x", "start_x", "end_x", "chr_y", "start_y", "end_y", "loopID", "Benign", "Malignant")
    ],
    file.path("Loops", "shared-loops.tsv"),
    sep = "\t",
    col.names = TRUE
)

# ==============================================================================
# Plots
# ==============================================================================
loginfo("Plotting figures")

# plots
gg_type <- (
    ggplot(data = collapsed_IDs_by_type, mapping = aes(x = Type))
    + geom_bar()
    + geom_text(stat = "count", mapping = aes(label=after_stat(count), vjust = -1))
    + scale_x_upset(name = "Tissue", order_by = "freq", n_intersections = Inf)
    + scale_y_continuous(
        limits = c(0, 5200),
        breaks = seq(0, 5000, by = 1000),
        name = "Loops"
    )
    + theme_minimal()
)
savefig(gg_type, file.path("Plots", "loop-calls.by-type"), width = 8, height = 12)

gg_sample <- (
    ggplot(data = collapsed_IDs_by_sample)
    + geom_bar(aes(x = SampleID))
    + labs(y = "Loops")
    + scale_x_upset(name = "Samples", order_by = "freq", n_intersections = 40)
    + theme_minimal()
)
savefig(gg_sample, file.path("Plots", "loop-calls.by-sample"), width = 40, height = 20)
