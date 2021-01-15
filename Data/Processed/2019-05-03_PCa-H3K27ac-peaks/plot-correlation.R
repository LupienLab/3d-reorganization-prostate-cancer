# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
suppressMessages(library("pheatmap"))

# ==============================================================================
# Data
# ==============================================================================
# load metadata
metadata <- fread("../../External/LowC_Samples_Data_Available.tsv", sep = "\t")
metadata <- metadata[Include == "Yes"]
SAMPLES <- paste0("PCa", metadata[, get("Sample ID")])
n <- length(SAMPLES)

batches <- fread("config.tsv", sep = "\t", header = TRUE)

# load long form correlation data
dt <- fread(
    "Peaks/correlation.tsv",
    header = FALSE,
    col.names = c("File_x", "File_y", "r"),
    sep = "\t"
)

# convert to wide form
mat <- as.data.frame(dcast(dt, File_x ~ File_y)[, 2:(n + 1)])
colnames(mat) <- SAMPLES
rownames(mat) <- SAMPLES

# ==============================================================================
# Analysis
# ==============================================================================
clust <- hclust(dist(mat, method = "euclidean"), method = "ward.D2")
clust_cut <- cutree(tree = clust, k = 2)
clust_cut <- data.table(
    Sample = SAMPLES,
    Patient = metadata[, get("Patient ID")],
    T2E = metadata[, get("T2E Status")],
    Cluster = clust_cut
)
clust_cut <- merge(
    x = clust_cut,
    y = batches[get("ChIP Type") == "H3K27ac", .SD, .SDcols = c("Patient ID", "Flowcell ID")],
    by.x = "Patient",
    by.y = "Patient ID"
)

print(clust_cut)

# ==============================================================================
# Plots
# ==============================================================================
