# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("DESeq2"))
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
    if (length(non_pos_cols) > 0) {
        mcols(gr) <- dt[, .SD, .SDcols = non_pos_cols]
    }
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

# load raw counts of H3K27ac pull down and input, then take differences
count_matrix <- NULL
library_sizes <- sapply(SAMPLES, function(s) 0)
for (s in SAMPLES) {
    dt_chip <- fread(
        file.path("Acetylation", paste0(s, "_H3K27ac.induced-region-counts.bed")),
        sep = "\t",
        header = FALSE,
        col.names = c("chr", "start", "end", "count", "supported", "width", "frac_supported")
    )
    dt_input <- fread(
        file.path("Acetylation", paste0(s, "_input.induced-region-counts.bed")),
        sep = "\t",
        header = FALSE,
        col.names = c("chr", "start", "end", "count", "supported", "width", "frac_supported")
    )
    # regions are in the same sorted order, so just taking difference
    dt_chip[, norm_count := pmax(0, count - dt_input$count)]
    count_matrix <- cbind(count_matrix, dt_chip$norm_count)
    library_sizes[s] <- sum(dt_chip$count)
}
colnames(count_matrix) <- SAMPLES

# load induced regions (the rows of the count matrix)
regions <- fread(
    "TAD-induced-regions.bed",
    sep = "\t",
    header = FALSE,
    col.names = c("chr", "start", "end")
)
regions_gr <- dt2gr(regions)


# ==============================================================================
# Analysis
# ==============================================================================
# create design matrices for each set of mut-vs-nonmut (this depends on the breakpoint being considered)
all_comparisons <- unique(tests$mutated_in)

# perform differential analysis for each test
for (mut_samples in all_comparisons[1]) {
    print(mut_samples)
    # create metadata table for samples
    meta <- data.table(Sample = SAMPLES, Mutated = "No", Size = library_sizes)
    meta[Sample %in% mut_samples, Mutated := "Yes"]
    meta[, Mutated := factor(Mutated, levels = c("No", "Yes"))]
    # create DESeqDataSet object
    dds <- DESeqDataSetFromMatrix(
        countData = count_matrix,
        colData = meta,
        design = ~ Mutated,
        rowRanges = regions_gr
    )
    # filter out regions with few reads in the region
    cat(">>filtering regions\n")
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep, ]
    cat(">>performing normalization\n")
    dds <- DESeq(dds)
    cat(">>extracting results\n")
    res <- results(dds)
    print(summary(res))
}


# ==============================================================================
# Plots
# ==============================================================================

