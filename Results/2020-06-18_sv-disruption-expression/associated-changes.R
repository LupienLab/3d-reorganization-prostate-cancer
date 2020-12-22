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

# combine all SVs that have some BND to compare inter-chromosomal events to strictly intra-chromosomal events
svtype_alteration_contingency_mat[, any_BND := BND + get("BND,DEL") + get("BND,DUP") + get("BND,UNKNOWN")]
svtype_alteration_contingency_mat[, only_same := DEL + DUP + INV + UNKNOWN]

svtype_alteration_contingency_mat <- svtype_alteration_contingency_mat[, .SD, .SDcols = c("altered_Exprs", "any_BND", "only_same")]

# perform goodness of fit tests
chisq.test(svtype_alteration_contingency_mat[, 2:3])

# contengency table for altered TAD and altered expression
tad_exprs_contingency <- sv_disruption[, .N, by = c("altered_TAD", "altered_Exprs")]
tad_exprs_contingency_mat <- dcast(
    tad_exprs_contingency,
    altered_TAD ~ altered_Exprs,
    value.var = "N",
    fill = 0
)

# perform independence test
chisq.test(tad_exprs_contingency_mat[, 2:3])


# ==============================================================================
# Plots
# ==============================================================================
# convert newly summed contingency table back to long form for plotting
svtype_alteration_contingency <- melt(
    svtype_alteration_contingency_mat,
    id.vars = "altered_Exprs",
    variable.name = "SV_Type",
    value.name = "N"
)
gg <- (
    ggplot(data = svtype_alteration_contingency)
    + geom_col(
        aes(x = SV_Type, y = N, fill = altered_Exprs),
        position = position_dodge()
    )
    + labs(x = "SV Type", y = "Count")
    + scale_x_discrete(
        breaks = c("only_same", "any_BND"),
        labels = c("Intra-chromosomal", "Inter-chromosomal")
    )
    + scale_fill_manual(
        name = "Altered Expression",
        breaks = c(FALSE, TRUE),
        labels = c("No", "Yes"),
        values = c("#AEA28E", "#FF6347")
    )
    + theme_minimal()
    + theme(
        legend.position = "bottom",
        panel.grid.minor.x = NULL,
        panel.grid.major.x = NULL
    )
)
ggsave(
    "Plots/sv-disruption-chromosomal.png",
    height = 12,
    width = 20,
    units = "cm"
)
ggsave(
    "Plots/sv-disruption-chromosomal.pdf",
    height = 12,
    width = 20,
    units = "cm"
)


# ==============================================================================
# Save data
# ==============================================================================
loginfo("Saving data")
fwrite(
    svtype_alteration_contingency_mat,
    "sv-disruption-tests.exprs.sv-type.tsv",
    sep = "\t"
)

fwrite(
    tad_exprs_contingency_mat,
    "sv-disruption-tests.exprs-TADs.tsv",
    sep = "\t"
)
