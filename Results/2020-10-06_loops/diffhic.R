# ==============================================================================
# Meta
# ==============================================================================
# diffHiC
# --------------------------------------
# Description: Perform differential analysis on the loop calls between tumours and benigns
# Author: James Hawley


# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("logging"))

loginfo("Loading packages")
suppressMessages(library("data.table"))
suppressMessages(library("diffHic"))
suppressMessages(library("regioneR"))
suppressMessages(library("Matrix"))

# ==============================================================================
# Functions
# ==============================================================================

# ==============================================================================
# Data
# ==============================================================================
loginfo("Loading data")
# load metadata
metadata <- fread("config.tsv", sep = "\t", header = TRUE)
metadata <- metadata[Include == "Yes", .SD]
SAMPLES <- list(
    "all" = metadata[, SampleID],
    "benign" = metadata[Type == "Benign", SampleID],
    "tumour" = metadata[Type == "Malignant", SampleID]
)

# load loop calls
loops <- fread("Loops/merged-loops.sample-counts.tsv")

# load bins and their IDs
bins <- toGRanges(
    fread("hg38.bins.bed", sep = "\t", header = FALSE, col.names = c("chr", "start", "end")),
    genome = "hg38"
)
bins <- trim(bins) # remove 0's from starting points in bins

sample_subset <- c(SAMPLES[["tumour"]][1:2], SAMPLES[["benign"]][1:2])

# load sparse matrices
matrices <- lapply(
    sample_subset,
    function(s) {
        fread(
            input = file.path("..", "2020-08-29_TADs-downsampled", "TMP", paste0(s, ".300000000.res_10000bp.chr1.mtx")),
            sep = "\t",
            header = TRUE
        )
    }    
)

names(matrices) <- sample_subset

# convert data tables to sparse matrices
sparse_matrices <- lapply(
    sample_subset,
    function(s) {
        sparseMatrix(
            i = matrices[[s]][, bin1_id + 1],
            j = matrices[[s]][, bin2_id + 1],
            x = matrices[[s]][, count],
            dims = rep(length(bins[seqnames(bins) == "chr1"]), 2)
        )
    }
)
names(sparse_matrices) <- sample_subset

# load into ContactMatrix objects
contact_matrices <- lapply(
    sample_subset,
    function(s) {
        ContactMatrix(
            matrix = sparse_matrices[[s]],
            anchor1 = bins[seqnames(bins) == "chr1", ],
            anchor2 = bins[seqnames(bins) == "chr1", ]
        )
    }
)
names(contact_matrices) <- sample_subset

# convert into InteractionSet object
data_combined <- mergeCMs(
    contact_matrices[[1]],
    contact_matrices[[2]],
    contact_matrices[[3]],
    contact_matrices[[4]],
    filter = 10
)
colnames(data_combined) <- sample_subset

# ==============================================================================
# Analysis
# ==============================================================================
loginfo("Performing calculations")

# list of all loop anchors for the ContactMatrix object
regions <- unique(c(
    toGRanges(
        loops[, .(chr = chr_x, start = start_x, end = end_x, anchor_ID = anchor_ID_x)],
        genome = "hg38"
    ),
    toGRanges(
        loops[, .(chr = chr_y, start = start_y, end = end_y, anchor_ID = anchor_ID_y)],
        genome = "hg38"
    )
))

# design matrix
design <- model.matrix(
    ~factor(metadata[SampleID %in% sample_subset, Type])
)
rownames(design) <- sample_subset
colnames(design) <- c("Intercept", "Tumour")

# only keep bins with > 10 reads across all samples
low_depth_idx <- which(rowSums(assay(data_combined)) >= 10)

# convert to DGEList object to use edgeR functions
y <- asDGEList(data_combined)

# ==============================================================================
# Save data
# ==============================================================================
loginfo("Saving data")