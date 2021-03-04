# ==============================================================================
# Meta
# ==============================================================================
# compartments
# ------------------------------------------------
# Author: James Hawley
# Description: Compartments analysis


# ==============================================================================
# Environment
# ==============================================================================
suppressWarnings(library("logging"))
loginfo("Loading packages")

suppressWarnings(library("data.table"))
suppressWarnings(library("ggplot2"))

CMPMT_DIR <- file.path(
    "..", "..", "Data", "Processed", "2019-06-18_PCa-LowC-sequencing", "Contacts"
)

# ==============================================================================
# Functions
# ==============================================================================


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
        dt[, pos := paste0(chrom, ":", start, "-", end)]
        dt[, bin := .I]
        dt[, Sample_ID := s]
        return(dt[, .SD, .SDcols = c("Sample_ID", "pos", "bin", "E1", "E2", "E3")])
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

# calculating GLM for compartments
model <- lm(E1 ~ as.factor(bin) + Type, data = eigs)

# alternative model formulation
mu_global <- eigs[, mean(E1, na.rm = TRUE)]
eig_pos <- eigs[, mean(E1 - mu_global, na.rm = TRUE), by = "pos"]
eig_pos <- eig_pos[complete.cases(eig_pos)]

# ==============================================================================
# Plots
# ==============================================================================
loginfo("Plotting figures")

gg <- (
    ggplot(data = eigs)
    + geom_density(
        mapping = aes(x = E1)
    )
    + facet_wrap(~ Sample_ID)
    + theme_minimal()
)
ggsave(
    "e1.png",
    gg,
    width = 12,
    height = 12,
    units = "cm"
)

# ==============================================================================
# Save Data
# ==============================================================================
loginfo("Saving data")

