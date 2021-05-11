# ==============================================================================
# Meta
# ==============================================================================
# cluster-by-compartment
# ------------------------------------------------
# Author: James Hawley
# Description: Cluster samples by compartment eigenvectors


# ==============================================================================
# Environment
# ==============================================================================
suppressWarnings(library("logging"))
loginfo("Loading packages")

suppressWarnings(library("data.table"))
suppressWarnings(library("ggplot2"))
suppressWarnings(library("pheatmap"))

CMPMT_DIR <- file.path(
    "..", "..", "Data", "Processed", "2019-06-18_PCa-LowC-sequencing", "Contacts"
)

# ==============================================================================
# Data
# ==============================================================================
loginfo("Loading data")

# load sample metadata
meta <- fread("config.tsv")
meta <- meta[Include == "Yes"]
SAMPLES <- meta[, Sample_ID]

# load compartment eigenvalues
eigs <- rbindlist(lapply(
    SAMPLES,
    function(s) {
        dt <- fread(
            file.path(CMPMT_DIR, paste0(s, ".compartments.cis.vecs.tsv")),
            sep = "\t",
            header = TRUE
        )
        dt[, bin_ID := .I]
        dt[, Sample_ID := s]
        return(dt[,
            .SD,
            .SDcols = c("Sample_ID", "chrom", "start", "end", "bin_ID", "E1", "E2", "E3")
        ])
    }
))

# add tissue type metadata to eigs
eigs <- merge(
    x = eigs,
    y = meta[, .SD, .SDcols = c("Sample_ID", "Label", "Type")],
    by = "Sample_ID"
)


# ==============================================================================
# Analysis
# ==============================================================================
loginfo("Calculating")

# convert to wide format
eigs_wide <- dcast(
    eigs,
    chrom + start + end ~ Sample_ID,
    value.var = "E1"
)
eigs_wide <- eigs_wide[complete.cases(eigs_wide)]

# convert to matrix for clustering
eigs_mtx <- as.matrix(eigs_wide[, .SD, .SDcols = 4:ncol(eigs_wide)])
colnames(eigs_mtx) <- colnames(eigs_wide)[4:ncol(eigs_wide)]
rownames(eigs_mtx) <- eigs_wide[, paste0(chrom, ":", start, "-", end)]

# calculating sample clustering
eigs_dist_mtx <- dist(t(eigs_mtx), method = "euclidean")

# ==============================================================================
# Plots
# ==============================================================================
loginfo("Plotting")
pheatmap(
    mat = eigs_dist_mtx,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    # show_colnames = TRUE,
    # show_rownames = TRUE,
    clustering_method = "ward.D2",
    filename = "compartments.heatmap.png",
    legend_breaks = meta[, Sample_ID],
    legend_labels = meta[, Label],
    labels_col = meta[, Label],
    labels_row = meta[, Label]
)
