# ==============================================================================
# Meta
# ==============================================================================
# T2E differential loops
# --------------------------------------
# Description: Assess the differential loops between T2E+ and T2E- primary tumours
# Author: James Hawley


# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("logging"))

loginfo("Loading packages")
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))

# ==============================================================================
# Functions
# ==============================================================================

# ==============================================================================
# Data
# ==============================================================================
loginfo("Loading data")
metadata <- fread("config.tsv", sep = "\t")
metadata <- metadata[Include == "Yes"]
SAMPLES <- list(
    "tumour" = metadata[Type == "Malignant", SampleID],
    "benign" = metadata[Type == "Benign", SampleID]
)

# read loop calls
loops <- fread("../2020-10-06_loops/Loops/merged-loops.tsv", sep = "\t")
loop_counts <- fread("../2020-10-06_loops/Loops/merged-loops.sample-counts.tsv", sep = "\t")

# ==============================================================================
# Analysis
# ==============================================================================
loginfo("Performing calculations")

# only keep high quality loops with a detection_scale > 2.5
hq_loop_IDs <- loop_counts[detection_scale > 2.5, loop_ID]
loops <- loops[loop_ID %in% hq_loop_IDs]

# get loop IDs by tumour type
loops <- merge(
    x = loops,
    y = metadata[, .SD, .SDcols = c("SampleID", "Type")],
    all = TRUE
)
typed_loops <- list(
    "benign" = loops[Type == "Benign", unique(loop_ID)],
    "tumour" = loops[Type == "Malignant", unique(loop_ID)]
)
typed_loops$intersection <- intersect(typed_loops$benign, typed_loops$tumour)
typed_loops$benign_exclusive <- setdiff(typed_loops$benign, typed_loops$intersection)
typed_loops$tumour_exclusive <- setdiff(typed_loops$tumour, typed_loops$intersection)
# count the "shared" loops as those that are called in any sample
# except if they are exclusive to one of the tissue types
typed_loops$shared <- setdiff(
    loops[, unique(loop_ID)],
    union(typed_loops$tumour_exclusive, typed_loops$benign_exclusive)
)

exclusive_loops <- list(
    "shared" = loop_counts[loop_ID %in% typed_loops$shared],
    "tumour" = loop_counts[loop_ID %in% typed_loops$tumour_exclusive],
    "benign" = loop_counts[loop_ID %in% typed_loops$benign_exclusive]
)

# ==============================================================================
# Save data
# ==============================================================================
loginfo("Saving data")

fwrite(
    exclusive_loops$tumour,
    "loops.tumour-specific.tsv",
    sep = "\t",
    col.names = TRUE
)
fwrite(
    exclusive_loops$benign,
    "loops.benign-specific.tsv",
    sep = "\t",
    col.names = TRUE
)
fwrite(
    exclusive_loops$shared,
    "loops.shared.tsv",
    sep = "\t",
    col.names = TRUE
)
