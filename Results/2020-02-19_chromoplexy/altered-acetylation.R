# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("DESeq2"))
suppressMessages(library("rsamtools"))
suppressMessages(library("GenomicRanges"))

# ==============================================================================
# Functions
# ==============================================================================
dt2gr <- function(dt) {
    '%ni%' <- Negate('%in%')
    non_pos_cols <- which(colnames(dt) %ni% c("chr", "start", "end"))
    gr <- GRanges(
        seqnames = dt$chr,
        ranges = IRanges(
            start = dt$start + 1,
            end = dt$end
        )
    )
    mcols(gr) <- dt[, .SD, .SDcols = non_pos_cols]
    return(gr)
}

split_comma_col <- function(v, f=identity) {
    splitv <- lapply(v, function(x) {strsplit(x, ",")[[1]]})
    return(lapply(splitv, f))
}

# ==============================================================================
# Data
# ==============================================================================
# load sample metadata
metadata <- fread(
    "../../Data/External/LowC_Samples_Data_Available.tsv",
    sep = "\t",
    header = TRUE
)
SAMPLES <- paste0("PCa", metadata[, get("Sample ID")])
metadata[, ChIP_file := paste0("../../Data/Processed/2019-05-03_PCa-H3K27ac-peaks/BAMs/Pca", get("Sample ID"), "_H3K27ac.sorted.dedup.bam")]
metadata[, Ctrl_file := paste0("../../Data/Processed/2019-05-03_PCa-H3K27ac-peaks/BAMs/Pca", get("Sample ID"), "_input.sorted.dedup.bam")]

# load breakpoints
breakpoints <- fread("sv-breakpoints.tsv", sep = "\t", header = TRUE)

# load TADs encompassing each breakpoint
tads <- fread("sv-breakpoints.TADs.tsv", sep = "\t", header = TRUE)
tads_gr <- dt2gr(tads)

# load tests metadata
tests <- fread("sv-disruption-tests.tsv", sep = "\t", header = TRUE)
# convert comma-separated values into lists
tests$breakpoint_indices <- split_comma_col(tests$breakpoint_indices, as.numeric)
tests$mutated_in <- split_comma_col(tests$mutated_in)

# create design matrices for each set of mut-vs-nonmut (this depends on the breakpoint being considered)
all_comparisons <- unique(tests$mutated_in])
all_designs <- lapply(
    all_comparisons,
    function(mut_samples) {
        dt <- data.table(Sample = SAMPLES, Mutated = FALSE)
        dt[Sample %in% mut_samples, Mutated := TRUE]
        design <- model.matrix(~ Mutated, dt)
        rownames(design) <- SAMPLES
        return(design)
    }
)

# ==============================================================================
# Analysis
# ==============================================================================

# ==============================================================================
# Plots
# ==============================================================================