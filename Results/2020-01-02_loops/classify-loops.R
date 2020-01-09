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
        "-p", "--prefix",
        type = "character",
        help = "Prefix for output files",
        default = "output"
    )
    ARGS <- PARSER$parse_args()
} else {
    ARGS = list(
        loops = file.path("cLoops_DefaultParameters", "PCa13266", "hic_loops_juicebox.txt"),
        cres = file.path("..", "..", "Data", "Processed", "2019-05-03_PCa-H3K27ac-peaks", "Peaks", "Pca13266_peaks.filtered.bedGraph"),
        prefix = "output"
    )
}

# ==============================================================================
# Functions
# ==============================================================================

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

# ==============================================================================
# Analysis
# ==============================================================================
# calculate intersection of loops and CREs
hits_x = as.data.table(findOverlaps(loops_x, peaks))
hits_y = as.data.table(findOverlaps(loops_y, peaks))

# calculate the number of loops where 0, 1, or 2 of its anchors overlap a CRE
intersecting_loop_hits = list(
    x_and_y = intersect(hits_x[, queryHits], hits_y[, queryHits])
)
intersecting_loop_hits[["x_and_noty"]] = setdiff(hits_x[, queryHits], intersecting_loop_hits[["x_and_y"]])
intersecting_loop_hits[["notx_and_y"]] = setdiff(hits_y[, queryHits], intersecting_loop_hits[["x_and_y"]])
intersecting_loop_hits[["notx_and_noty"]] = setdiff(
    1:loops_dt[, .N],
    unlist(intersecting_loop_hits)
)
# map list elements back to the loop IDs and store results in `loops_dt`
intersecting_loop_IDs = lapply(
    names(intersecting_loop_hits),
    function(intersecting_class) {
        loops_dt[intersecting_loop_hits[[intersecting_class]], V1 := intersecting_class]
        return(loops_x[intersecting_loop_hits[[intersecting_class]]]$LoopID)
    }
)
names(intersecting_loop_IDs) = names(intersecting_loop_hits)
colnames(loops_dt)[16] = "CRE_Overlap"

# calculate the number of CREs overlapping at least 1 anchor in a loop pair
intersecting_CRE_hits = list(
    in_loop = union(hits_x[, subjectHits], hits_y[, subjectHits])
)
intersecting_CRE_hits$no_loop = setdiff(1:peaks_dt[, .N], intersecting_CRE_hits$in_loop)
# sanity check: ensure that a loop isn't called on a single element
if (length(unique(intersecting_CRE_hits$in_loop)) != length(intersecting_CRE_hits$in_loop)) {
    print("Loops called on a single element in this sample")
}

# map list elements back to the CRE IDs and store results in `peaks_dt`
intersecting_CRE_IDs = lapply(
    names(intersecting_CRE_hits),
    function(intersecting_class) {
        peaks_dt[intersecting_CRE_hits[[intersecting_class]], V1 := intersecting_class]
        return(peaks[intersecting_CRE_hits[[intersecting_class]]]$PeakID)
    }
)
names(intersecting_CRE_IDs) = names(intersecting_CRE_hits)
colnames(peaks_dt)[6] = "Loop_Overlap"
peaks_dt[Loop_Overlap == "in_loop", Loop_Overlap := TRUE]
peaks_dt[Loop_Overlap == "no_loop", Loop_Overlap := FALSE]

# ==============================================================================
# Save data
# ==============================================================================
# remove "colour" column and save to TSV
loops_dt[, colour := NULL]
fwrite(loops_dt, paste(ARGS$prefix, "loops", "tsv", sep = "."), sep = "\t", col.names = TRUE)
fwrite(peaks_dt, paste(ARGS$prefix, "peaks", "tsv", sep = "."), sep = "\t", col.names = TRUE)
