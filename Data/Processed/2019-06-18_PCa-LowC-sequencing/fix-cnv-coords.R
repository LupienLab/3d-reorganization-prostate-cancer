# ==============================================================================
# Environment
# ==============================================================================
suppressWarnings(library("argparse"))

if (!interactive()) {
	parser <- ArgumentParser(description = "Fix coordinates from HiNT output segmentations")
	parser$add_argument(
		"bed",
		type = "character",
		help = "BED-like file"
	)
	parser$add_argument(
		"sizes",
		type = "character",
		help = "Genome chromosome sizes file"
	)
	parser$add_argument(
		"-o", "--output",
		type = "character",
		help = "Output file"
	)
	args <- parser$parse_args()
} else {
	args <- list(
		"bed" = file.path("CNV", "segmentation", "PCa13266_CNV_segments.txt"),
		"sizes" = "hg38.sizes.txt",
		output = file.path("CNV", "segmentation", "PCa13266_CNV_segments.corrected.txt")
	)
}

suppressWarnings(library("data.table"))

# ==============================================================================
# Functions
# ==============================================================================
correct_pos <- function(dt) {
	# merge genome chromosome sizes into table
	dt <- merge(
		x = dt,
		y = sizes,
		by = "chr",
		all.x = TRUE
	)
	# cap position values as required
	dt[, start := as.integer(pmax(0, start))]
	dt[, end := as.integer(pmax(0, pmin(end, chrom_size)))]
	return(dt[, .SD, .SDcols = -"chrom_size"])
}


# ==============================================================================
# Data
# ==============================================================================
# load chromosome sizes
sizes <- fread(args$sizes, sep = "\t", col.names = c("chr", "chrom_size"))

# load file to correct
segments <- fread(
	args$bed,
	sep = "\t",
	col.names = c("chr", "start", "end", "binNum", "observed", "expected", "log2CopyRatio", "p")
)


# ==============================================================================
# Analysis
# ==============================================================================
# fix off-by-1 position for end coordinate
segments[, end := end + 1]
# correct positions as needed
segments_corrected <- correct_pos(segments)

# ==============================================================================
# Save data
# ==============================================================================
fwrite(
	segments_corrected,
	args$output,
	sep = "\t",
	col.names = TRUE
)

