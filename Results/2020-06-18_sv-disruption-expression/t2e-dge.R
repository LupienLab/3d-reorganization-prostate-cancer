# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("rlang"))
suppressMessages(library("vctrs"))
suppressMessages(library("sleuth"))
suppressMessages(library("data.table"))

KALLISTO_DIR <- file.path("..", "..", "Data", "Processed", "2020-06-17_PCa-RNA-seq")

# ==============================================================================
# Functions
# ==============================================================================

# ==============================================================================
# Data
# ==============================================================================
# load sample metadata
metadata <- fread(
    "../../Data/External/LowC_Samples_Data_Available.tsv",
    sep = "\t",
    header = TRUE
)
metadata <- metadata[Include == "Yes"]
metadata[, SampleID :=  paste0("PCa", get("Sample ID"))]
SAMPLES <- metadata$SampleID

# add kallisto paths
metadata[, Path := file.path(KALLISTO_DIR, SampleID)]

design <- metadata[, .(sample = SampleID, condition = get("T2E Status"), path = Path)]


# ==============================================================================
# Analysis
# ==============================================================================
# create sleuth object
so <- sleuth_prep(design, extra_bootstrap_summary = TRUE)

# fit full model
so <- sleuth_fit(so, ~condition, "full")

# fit reduced (null) model
so <- sleuth_fit(so, ~1, "reduced")

# perform differential analysis
so <- sleuth_lrt(so, "reduced", "full")

# save data
saveRDS(so, "sleuth-object.rds")

# extract results
sleuth_table <- as.data.table(sleuth_results(so, "reduced:full", "lrt", show_all = FALSE))
fwrite(sleuth_table, "t2e-comparison.tsv", sep = "\t", col.names = TRUE)
