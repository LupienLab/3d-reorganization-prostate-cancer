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
    "t2e" = metadata[T2E == "Yes", SampleID],
    "nont2e" = metadata[T2E == "No", SampleID]
)

# read loop calls
loops <- fread("../2020-10-06_loops/Loops/merged-loops.tsv", sep = "\t")
loop_counts <- fread("../2020-10-06_loops/Loops/merged-loops.sample-counts.tsv", sep = "\t")

# ==============================================================================
# Analysis
# ==============================================================================
loginfo("Performing calculations")
# start with the loops that are called in at least 3 samples
consistent_loops <- loop_counts[Malignant >= 3, loopID]
loops <- loops[loopID %in% consistent_loops]

# get loop IDs by tumour type
loops <- merge(loops, metadata[, .SD, .SDcols = c("SampleID", "T2E")])
typed_loops <- list(
    "t2e" = loops[T2E == "Yes", unique(loopID)],
    "nont2e" = loops[T2E == "No", unique(loopID)]
)
typed_loops$intersection <- intersect(typed_loops$t2e, typed_loops$nont2e)
typed_loops$nont2e_exclusive <- setdiff(typed_loops$nont2e, typed_loops$intersection)
typed_loops$t2e_exclusive <- setdiff(typed_loops$t2e, typed_loops$intersection)

exclusive_loops <- list(
    "shared" = loop_counts[loopID %in% typed_loops$intersection],
    "t2e" = loop_counts[loopID %in% typed_loops$t2e_exclusive],
    "nont2e" = loop_counts[loopID %in% typed_loops$nont2e_exclusive]
)

# ==============================================================================
# Save data
# ==============================================================================
loginfo("Saving data")

fwrite(
    exclusive_loops$t2e,
    "loops.T2E-specific.tsv",
    sep = "\t",
    col.names = TRUE
)
fwrite(
    exclusive_loops$nont2e,
    "loops.nonT2E-specific.tsv",
    sep = "\t",
    col.names = TRUE
)
fwrite(
    exclusive_loops$shared,
    "loops.shared.tsv",
    sep = "\t",
    col.names = TRUE
)
