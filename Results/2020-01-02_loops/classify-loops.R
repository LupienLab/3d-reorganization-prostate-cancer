# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("GenomicRanges"))
suppressMessages(library("argparse"))

if (!interactive()) {
    PARSER <- argparse::ArgumentParser(
        description = "Classify called loops according to their location"
    )
    PARSER$add_argument(
        "loops",
        type = "character",
        help = "BEDPE file of called loops"
    )
    PARSER$add_argument(
        "cres",
        type = "character",
        help = "BED file containing cis-regulatory element loci"
    )
    PARSER$add_argument(
        "boundaries",
        type = "character",
        help = "BED file containing TAD boundaries across window sizes"
    )
    PARSER$add_argument(
        "-p", "--prefix",
        type = "character",
        help = "Prefix for output files",
        default = "output"
    )
    ARGS <- PARSER$parse_args()
} else {
    ARGS = list(
        loops = file.path("Loops", "PCa13266", "hic_loops_juicebox.txt"),
        cres = file.path("..", "..", "Data", "Processed", "2019-05-03_PCa-H3K27ac-peaks", "Peaks", "Pca13266_peaks.filtered.bedGraph"),
        boundaries = file.path("..", "2020-01-15_TAD-aggregation", "resolved-TADs", "PCa13266.40000bp.aggregated-boundaries.tsv"),
        prefix = file.path("Classification", "PCa13266")
    )
}

# ==============================================================================
# Data
# ==============================================================================
# load loop calls
loops_dt = fread(
    ARGS$loops,
    sep = "\t",
    header = TRUE,
    col.names = c(
        "chr_x", "start_x", "end_x",
        "chr_y", "start_y", "end_y",
        "colour", "observed", "LoopID", "FDR", "EnrichmentScore", "Distance",
        "log10p_binom", "log10p_poisson", "log10p_hypg"
    )
)

# convert to GRanges, splitting into pairs
loops_x = GRanges(
    seqnames = loops_dt$chr_x,
    ranges = IRanges(
        start = loops_dt$start_x + 1,
        end = loops_dt$end_x
    ),
    LoopID = loops_dt$LoopID
)
loops_y = GRanges(
    seqnames = loops_dt$chr_y,
    ranges = IRanges(
        start = loops_dt$start_y + 1,
        end = loops_dt$end_y
    ),
    LoopID = loops_dt$LoopID
)

# load cis-regulatory element locations (defined by H3K27ac peaks)
peaks_dt = fread(ARGS$cres, sep = "\t", header = FALSE, col.names = c("chr", "start", "end", "tenlog10q"))

# add unique ID for each peak to keep track
peaks_dt[, PeakID := paste0(chr, ":", start, "-", end)]

# convert to GRanges
peaks = GRanges(
    seqnames = peaks_dt$chr,
    ranges = IRanges(
        start = peaks_dt$start + 1,
        end = peaks_dt$end
    ),
    PeakID = peaks_dt$PeakID
)

# load TAD boudnaries
boundaries_dt = fread(ARGS$boundaries, sep = "\t", header = TRUE)
boundaries_dt[, BoundID := paste0(chr, ":", pos)]
boundaries = GRanges(
    seqnames = boundaries_dt$chr,
    ranges = IRanges(start = boundaries_dt$pos + 1 - 20000, width = 40000),
    BoundID = boundaries_dt$BoundID
)

# ==============================================================================
# Analysis
# ==============================================================================
# Calculate intersection of loops and CREs
# --------------------------------------
peak_hits_x = as.data.table(findOverlaps(loops_x, peaks))
peak_hits_y = as.data.table(findOverlaps(loops_y, peaks))

# calculate the number of loops where 0, 1, or 2 of its anchors overlap a CRE
loops_dt[, Overlapping_Loop_x := FALSE]
loops_dt[, Overlapping_Loop_y := FALSE]
loops_dt[peak_hits_x[, queryHits], Overlapping_Loop_x := TRUE]
loops_dt[peak_hits_y[, queryHits], Overlapping_Loop_y := TRUE]
loops_dt[, CRE_Overlap := Overlapping_Loop_x + Overlapping_Loop_y]

# add Peak IDs for loops where there are overlapping peaks
x_loop_peakIDs = peaks_dt[peak_hits_x[, subjectHits], PeakID]
y_loop_peakIDs = peaks_dt[peak_hits_y[, subjectHits], PeakID]
loops_dt[peak_hits_x[, queryHits], Overlapping_Loop_x_ID:= x_loop_peakIDs]
loops_dt[peak_hits_y[, queryHits], Overlapping_Loop_y_ID:= y_loop_peakIDs]

# calculate the number of CREs overlapping at least 1 anchor in a loop pair
idx_in_loop = unique(union(peak_hits_x[, subjectHits], peak_hits_y[, subjectHits]))
peaks_dt[, Loop_Overlap := FALSE]
peaks_dt[idx_in_loop, Loop_Overlap := TRUE]

# Calculate intersection of loops and TAD boundaries
# --------------------------------------
bound_hits_x = as.data.table(findOverlaps(loops_x, boundaries))
bound_hits_y = as.data.table(findOverlaps(loops_y, boundaries))

# calculate the number of loops where 0, 1, or 2 of its anchors overlap a CRE
loops_dt[, Overlapping_Boundary_x := FALSE]
loops_dt[, Overlapping_Boundary_y := FALSE]
loops_dt[bound_hits_x[, queryHits], Overlapping_Boundary_x := TRUE]
loops_dt[bound_hits_y[, queryHits], Overlapping_Boundary_y := TRUE]
loops_dt[, TAD_Overlap := Overlapping_Boundary_x + Overlapping_Boundary_y]

# add Peak IDs for loops where there are overlapping peaks
x_loop_boundIDs = boundaries_dt[bound_hits_x[, subjectHits], BoundID]
y_loop_boundIDs = boundaries_dt[bound_hits_y[, subjectHits], BoundID]
boundaries_dt[bound_hits_x[, queryHits], Overlapping_Loop_x_ID:= x_loop_boundIDs]
boundaries_dt[bound_hits_y[, queryHits], Overlapping_Loop_y_ID:= y_loop_boundIDs]

# calculate the number of CREs overlapping at least 1 anchor in a loop pair
idx_in_loop = unique(union(bound_hits_x[, subjectHits], bound_hits_y[, subjectHits]))
boundaries_dt[, Loop_Overlap := FALSE]
boundaries_dt[idx_in_loop, Loop_Overlap := TRUE]

# ==============================================================================
# Save data
# ==============================================================================
# remove "colour" column and save to TSV
loops_dt[, colour := NULL]
fwrite(loops_dt, paste(ARGS$prefix, "loops", "tsv", sep = "."), sep = "\t", col.names = TRUE)
fwrite(loops_dt, paste(ARGS$prefix, "bounds", "tsv", sep = "."), sep = "\t", col.names = TRUE)
fwrite(peaks_dt, paste(ARGS$prefix, "peaks", "tsv", sep = "."), sep = "\t", col.names = TRUE)
