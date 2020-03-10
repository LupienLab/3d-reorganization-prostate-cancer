# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("DiffBind"))

# ==============================================================================
# Functions
# ==============================================================================

# ==============================================================================
# Data
# ==============================================================================
# load metadata
metadata <- fread("../../External/LowC_Samples_Data_Available.tsv", sep = "\t")
SAMPLES = paste0("PCa", metadata[, get("Sample ID")])
n = length(SAMPLES)

batches <- fread("config.tsv", sep = "\t", header = TRUE)

# merge flowcell IDs
metadata <- merge(
    x = metadata,
    y = batches[get("ChIP Type") == "H3K27ac", .SD, .SDcols = c("Patient ID", "Flowcell ID")],
    by = "Patient ID",
)
colnames(metadata)[11] <- "ChIP_Flowcell"
metadata <- merge(
    x = metadata,
    y = batches[get("ChIP Type") == "input", .SD, .SDcols = c("Patient ID", "Flowcell ID")],
    by = "Patient ID",
)
colnames(metadata)[12] <- "Input_Flowcell"

mod <- metadata[, model.matrix(get("Patient ID") ~ ChIP_Flowcell + Input_Flowcell)]
# ==============================================================================
# Analysis
# ==============================================================================
# configuration table for condition comparisons
config <- data.table(
    SampleID = metadata[, get("Patient ID")],
    Tissue = "Prostate",
    Factor = metadata[, ChIP_Flowcell],
    Condition = metadata[, get("T2E Status")],
    bamReads = Sys.glob(file.path("..", "..", "Raw", paste0("*_", metadata[, ChIP_Flowcell]), "Aligned", metadata[, paste0("Pca", get("Sample ID"), "_H3K27ac.filtered.dedup.sorted.bam")])),
    controlReads = Sys.glob(file.path("..", "..", "Raw", paste0("*_", metadata[, Input_Flowcell]), "Aligned", metadata[, paste0("Pca", get("Sample ID"), "_input.filtered.dedup.sorted.bam")]))
)

db <- dba()

# ==============================================================================
# Plots
# ==============================================================================