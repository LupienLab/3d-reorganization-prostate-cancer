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

MIN_SAMPLES <- 2

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
# add loop_ID as a metadata column
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
    loop_ID = paste0(subjectHits(hits_x), "_", subjectHits(hits_y))
)]

# merge in tissue type information into the loops
loops <- merge(
    x = loops,
    y = metadata[, .SD, .SDcols = c("SampleID", "Type")],
    by = "SampleID"
)

# some loops, identified separately in the same sample, may now be merged because of loop calls from other samples
# this step ensures that the single loop isn't double-counted
loops <- loops[,
    .(
        start_x = min(start_x),
        end_x = max(end_x),
        start_y = min(start_y),
        end_y = max(end_y),
        fdr = min(fdr),
        detection_scale = max(detection_scale)
    ),
    by = c(
        "SampleID", "chr_x", "chr_y",
        "anchorID_x", "anchorID_y", "loop_ID", "Type"
    )
]

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
    "loop_ID" = loops[, loop_ID],
    "detection_scale" = loops[, detection_scale],
    "fdr" = loops[, fdr]
)

unique_loops <- merged_loops[,
    .(
        chr_x = unique(chr_x),
        start_x = unique(start_x),
        end_x = unique(end_x),
        chr_y = unique(chr_y),
        start_y = unique(start_y),
        end_y = unique(end_y),
        anchor_ID_x = unique(anchor_ID_x),
        anchor_ID_y = unique(anchor_ID_y),
        fdr = min(fdr),
        detection_scale = max(detection_scale)
    ),
    keyby = "loop_ID"
]

# count the number of times a loop occurs
loop_counts <- dcast(
    # count loops byt their occurrences in each tissue type
    merged_loops[, .N, keyby = c("loop_ID", "Type")],
    # pivot the table to a wide format, rows are loops, columns are counts in benign and tumour
    loop_ID ~ Type,
    value.var = "N",
    # fill any NAs as 0
    fill = 0
)

# merge back in position information for loops
unique_loops <- merge(
    x = unique_loops,
    y = loop_counts,
    by = "loop_ID",
    all = TRUE
)


# ==============================================================================
# Save data
# ==============================================================================
loginfo("Saving data")
fwrite(
    merged_loops,
    file.path("Loops", "merged-loops.tsv"),
    sep = "\t",
    col.names = TRUE
)
fwrite(
    # put the columns in an order compatible with BEDPE format
    unique_loops,
    file.path("Loops", "merged-loops.sample-counts.tsv"),
    sep = "\t",
    col.names = TRUE
)

# save benign-specific loops
fwrite(
    # put the columns in an order compatible with BEDPE format
    unique_loops[
        Benign >= MIN_SAMPLES & Malignant == 0,
        .SD,
        .SDcols = c(
            "chr_x", "start_x", "end_x",
            "chr_y", "start_y", "end_y",
            "anchor_ID_x", "anchor_ID_y", "loop_ID",
            "fdr", "detection_scale"
        )
    ],
    file.path("Loops", "benign-specific-loops.tsv"),
    sep = "\t",
    col.names = TRUE
)

# save tumour-specific loops
fwrite(
    # put the columns in an order compatible with BEDPE format
    unique_loops[
        Benign == 0 & Malignant >= MIN_SAMPLES,
        .SD,
        .SDcols = c(
            "chr_x", "start_x", "end_x",
            "chr_y", "start_y", "end_y",
            "anchor_ID_x", "anchor_ID_y", "loop_ID",
            "fdr", "detection_scale"
        )
    ],
    file.path("Loops", "tumour-specific-loops.tsv"),
    sep = "\t",
    col.names = TRUE
)

# save shared tumour-benign loops
fwrite(
    # put the columns in an order compatible with BEDPE format
    unique_loops[
        Benign > 0 & Malignant > 0,
        .SD,
        .SDcols = c(
            "chr_x", "start_x", "end_x",
            "chr_y", "start_y", "end_y",
            "anchor_ID_x", "anchor_ID_y", "loop_ID",
            "fdr", "detection_scale",
            "Benign", "Malignant"
        )
    ],
    file.path("Loops", "shared-loops.tsv"),
    sep = "\t",
    col.names = TRUE
)
