# ==============================================================================
# Meta
# ==============================================================================
# Get loop anchors and unify them across a given sample
# --------------------------------------
# Description: Take 2D loop calls and identify the unique anchors
# Author: James Hawley


# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("logging"))
suppressMessages(library("argparse"))

if (!interactive()) {
	PARSER <- argparse::ArgumentParser(
		description = "Take 2D loop calls and identify the unique anchors"
	)
	PARSER$add_argument(
		"prefix",
		type = "character",
		help = "File prefix, used for both input and output"
	)
	ARGS <- PARSER$parse_args()
} else {
	ARGS <- list(
		"prefix" = file.path("..", "..", "results", "Loops", "PCa13266"),
	)
}


loginfo("Loading packages")
suppressMessages(library("data.table"))
suppressMessages(library("GenomicRanges"))
suppressMessages(library("regioneR"))
suppressMessages(library("ggplot2"))
suppressMessages(library("ggupset"))

CHRS <- paste0("chr", c(1:22, "X", "Y"))
DEPTHS <- 50000000 * (1:6)
REPS <- 1:10

RES_DIR <- file.path("..", "..", "results", "2021-07-19_loops-downsampled")
LOOP_DIR <- file.path(RES_DIR, "Loops")


# ==============================================================================
# Data
# ==============================================================================
loginfo("Loading data")

# load loop calls
loops <- rbindlist(lapply(
	DEPTHS,
	function(depth) {
		dt1 <- rbindlist(lapply(
			REPS,
			function(r) {
				dt2 <- fread(
					paste0(
						ARGS$prefix,
						".depth_", as.integer(depth),
						".res_10000bp.",
						"rep_", r,
						".loops.tsv"
					),
					sep = "\t",
					header = TRUE
				)
				dt2[, Replicate := r]
				return(dt2)
			}
		))
		dt1[, Seq_Depth := depth]
		return(dt1)
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

# some loops, identified separately in the same depth, may now be merged because of loop calls from other depths
# this step ensures that the single loop isn't double-counted
loops <- loops[,
	.(
		start_x = min(start_x),
		end_x = max(end_x),
		start_y = min(start_y),
		end_y = max(end_y),
		fdr = min(fdr),
		detection_scale = max(detection_scale),
	),
	by = c(
		"chr_x", "chr_y",
		"Replicate", "Seq_Depth",
		"anchorID_x", "anchorID_y", "loop_ID"
	)
]

# merge resulting loop calls to give a set of loops based on the anchors from all samples
corrected_anchors_x <- loop_anchors[loops[, anchorID_x]]
corrected_anchors_y <- loop_anchors[loops[, anchorID_y]]

# combine all loop calls together into a single table
merged_loops <- data.table(
	"chr_x" = as.vector(seqnames(corrected_anchors_x)),
	"start_x" = start(corrected_anchors_x),
	"end_x" = end(corrected_anchors_x),
	"chr_y" = as.vector(seqnames(corrected_anchors_y)),
	"start_y" = start(corrected_anchors_y),
	"end_y" = end(corrected_anchors_y),
	"Replicate" = loops[, as.integer(Replicate)],
	"Seq_Depth" = loops[, as.integer(Seq_Depth)],
	"anchor_ID_x" = loops[, as.integer(anchorID_x)],
	"anchor_ID_y" = loops[, as.integer(anchorID_y)],
	"loop_ID" = loops[, loop_ID],
	"fdr" = loops[, fdr],
	"detection_scale" = loops[, detection_scale]
)


# ==============================================================================
# Save data
# ==============================================================================
loginfo("Saving data")
fwrite(
	merged_loops,
	paste0(ARGS$prefix, ".all-iterations.loops.tsv"),
	sep = "\t",
	col.names = TRUE
)

