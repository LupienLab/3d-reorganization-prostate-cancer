# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("scales"))
suppressMessages(library("ggplot2"))
suppressMessages(library("Matrix"))

CHRS <- paste0("chr", c(1:22, "X"))

# ==============================================================================
# Functions
# ==============================================================================
#' Save figures in multiple formats
#'
#' @param gg ggplot object
#' @param prefix Prefix for output file
#' @param ext Output extensions
#' @param dpi DPI resolution
savefig = function(gg, prefix, ext = c("png", "pdf"), width = 20, height = 12, dpi = 400) {
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

