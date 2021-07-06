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

scc_mat <- as.matrix(scc[, .SD, .SDcols = -"SampleID"])
rownames(scc_mat) <- SAMPLES
colnames(scc_mat) <- SAMPLES

# ==============================================================================
# Plots
# ==============================================================================
loginfo("Plotting heatmap")

pheatmap(
	mat = scc_mat,
	cluster_rows = TRUE,
	cluster_cols = TRUE,
	filename = file.path(DIR[["plot"]], "scc.heatmap.png")
)


# ==============================================================================
# Save data
# ==============================================================================
loginfo("Saving data")
fwrite(
	agg_data,
	file.path(DIR[["res"]], "scc.tsv"),
	sep = "\t",
	col.names = TRUE
)


