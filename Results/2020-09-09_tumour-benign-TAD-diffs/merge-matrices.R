# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("scales"))
suppressMessages(library("ggplot2"))
suppressMessages(library("Matrix"))

CHRS <- paste0("chr", c(1:22, "X"))


# ==============================================================================
# Data
# ==============================================================================
cat("Loading Data\n")
# load metadata
metadata <- fread("config.tsv", sep = "\t", header = TRUE)
metadata <- metadata[Include == "Yes"]
TUMOUR_SAMPLES <- metadata[Type == "Malignant", SampleID]
BENIGN_SAMPLES <- metadata[Type == "Benign", SampleID]


# ==============================================================================
# Analysis
# ==============================================================================
# loop over each chromosome
for (chrom in rev(CHRS)) {
    tumour_mats <- lapply(
        TUMOUR_SAMPLES,
        function(id) {
            as.matrix(fread(
                cmd = paste0("zcat TMP/", id, ".300000000.res_40000bp.", chrom, ".noheader.mtx.gz"),
                sep = "\t"
            ))
        }
    )
    names(tumour_mats) <- TUMOUR_SAMPLES
    for (s in tumour_mats) {print(head(tumour_mats[[s]]))}
    readline()
}


# ==============================================================================
# Save data
# ==============================================================================

