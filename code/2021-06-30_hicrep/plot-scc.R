# ==============================================================================
# Meta
# ==============================================================================
# Name
# --------------------------------------
# Description: Plot structural similarity coefficient values
# Author: James Hawley


# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("logging"))

loginfo("Loading packages")
suppressMessages(library("data.table"))
suppressMessages(library("Matrix"))
suppressMessages(library("pheatmap"))

DIR <- list(
	"res" = file.path("..", "..", "results", "2021-06-30_hicrep")
)
DIR[["plot"]] <- file.path(DIR[["res"]], "Plots")


# ==============================================================================
# Data
# ==============================================================================
loginfo("Loading data")

# metadata
meta <- fread("config.tsv", sep = "\t", header = TRUE)
meta <- meta[(Include == "Yes") & (Source == "Primary")]
SAMPLES <- meta[, SampleID]

# SCC values
scc <- fread(
	file.path(DIR[["res"]], "scc.tsv"),
	sep = "\t",
	header = TRUE
)

# convert the values into a matrix for plotting
scc_mat <- as.matrix(scc[, .SD, .SDcols = -"SampleID"])
rownames(scc_mat) <- SAMPLES
colnames(scc_mat) <- SAMPLES

# the SCC values will not be exactly equal because of the downsampling
# to make sure the clustering is the same for the rows and columns, we
# symmetrize the matrix. The true differences between the two are
# negligible, compared to the differences between samples
# (i.e. ||skew_mat|| < 1e-4)
symm_mat <- symmpart(scc_mat)
skew_mat <- skewpart(scc_mat)


# ==============================================================================
# Plots
# ==============================================================================
loginfo("Plotting heatmap")

pheatmap(
	mat = symm_mat,
	cluster_rows = TRUE,
	cluster_cols = TRUE,
	labels_row = meta[, Label],
	labels_col = meta[, Label],
	filename = file.path(DIR[["plot"]], "scc.heatmap.png")
)
pheatmap(
	mat = symm_mat,
	cluster_rows = TRUE,
	cluster_cols = TRUE,
	labels_row = meta[, Label],
	labels_col = meta[, Label],
	filename = file.path(DIR[["plot"]], "scc.heatmap.pdf")
)
