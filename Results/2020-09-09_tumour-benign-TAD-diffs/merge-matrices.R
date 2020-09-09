# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("MASS"))


CHRS <- paste0("chr", c(1:22, "X"))


# ==============================================================================
# Functions
# ==============================================================================
# simplify the process of reading in exported matrix files
read_contacts <- function(id, chrom) {
    # read zipped matrix file
    dt <- fread(
        cmd = paste0("zcat TMP/", id, ".300000000.res_40000bp.", chrom, ".noheader.mtx.gz"),
        sep = "\t"
    )
    # for some reason the last column is filled with NAs, but these shouldn't be here
    na_cols <- dt[, which(apply(.SD, 2, function(x) all(is.na(x))))]
    # remove these columns
    dt <- dt[, .SD, .SDcols = -na_cols]
    # convert to a matrix and remove column names
    mat <- as.matrix(dt)
    colnames(mat) <- NULL
    # return the matrix itself
    return(mat)
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
chrom_mats <- list(
    "tumour" = list(),
    "benign" = list()
)

for (chrom in CHRS) {
    cat(chrom, "\n")
    tumour_mats <- lapply(TUMOUR_SAMPLES[1:2], read_contacts, chrom)
    # calculate the mean contact matrix from all the tumour matrices
    chrom_mats[["tumour"]][[chrom]] <- Reduce("+", tumour_mats) / length(tumour_mats)
    benign_mats <- lapply(BENIGN_SAMPLES[1:2], read_contacts, chrom)
    # calculate the mean contact matrix from all the tumour matrices
    chrom_mats[["benign"]][[chrom]] <- Reduce("+", benign_mats) / length(benign_mats)
}


# ==============================================================================
# Save data
# ==============================================================================
for (chrom in CHRS) {
    cat(chrom, "\n")
    for (i in c("tumour", "benign")) {
        cat("\t", i, "\n")
        dt <- as.data.table(chrom_mats[[i]][[chrom]])
        fwrite(
            dt,
            paste0("TMP/", i, ".300000000.res_40000bp.", chrom, ".mtx"),
            sep = "\t",
            col.names = FALSE
        )
    }
}

