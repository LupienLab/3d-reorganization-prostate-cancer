# ==============================================================================
# Meta
# ==============================================================================
# associated-changes
# --------------------------------------
# Description: Calculate associated changes with SVs and expression
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

# load SVs and impact on expression
sv_disruption <- fread("summary-sv-disruption-impact.tsv")
# ==============================================================================
# Analysis
# ==============================================================================
loginfo("Performing calculations")

# calculate contingency table for altered status based on SV breakpoint type
svtype_alteration_contingency <- sv_disruption[, .N, by = c("altered_Exprs", "SV_type")]
svtype_alteration_contingency_mat <- dcast(
    svtype_alteration_contingency,
    altered_Exprs ~ SV_type,
    value.var = "N",
    fill = 0
)

# perform goodness of fit tests
chisq.test(svtype_alteration_contingency_mat[, 2:9])

# ==============================================================================
# Save data
# ==============================================================================
loginfo("Saving data")
fwrite(
    svtype_alteration_contingency_mat,
    "sv-disruption-tests.exprs.sv-type.tsv",
    sep = "\t"
)
