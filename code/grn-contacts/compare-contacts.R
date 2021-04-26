# ==============================================================================
# Meta
# ==============================================================================
# compare-contacts
# ------------------------------------------------
# Author: James Hawley
# Description: Compare obs/exp contact frequencies from GRNs


# ==============================================================================
# Environment
# ==============================================================================
suppressWarnings(library("logging"))
loginfo("Loading packages")

suppressWarnings(library("data.table"))
suppressWarnings(library("ggplot2"))

RESULT_DIR = file.path("..", "..", "Results", "grn-contacts")
EXPRS_DIR = file.path("..", "..", "Results", "2020-06-18_sv-disruption-expression")
PLOT_DIR = file.path(RESULT_DIR, "Plots")


# ==============================================================================
# Data
# ==============================================================================
loginfo("Loading data")

# load contact frequencies
contacts <- fread(
    file.path(RESULT_DIR, "contacts-extracted.tsv"),
    sep = "\t",
    header = TRUE,
    col.names = c(
        "contact_Sample_ID",
        "contact_event_ID",
        "query_gene_id",
        "query_gene_name",
        "query_prom_chr",
        "query_prom_start",
        "query_prom_end",
        "enh_target_gene_id",
        "enh_chr",
        "enh_start",
        "enh_end",
        "sample_mutated_in_event",
        "Mean_Obs_Exp_Freq"
    )
)

# load differential expression results
exprs <- fread(
    file.path(EXPRS_DIR, "summary-sv-disruption.tsv"),
    sep = "\t",
    header = TRUE
)


# ==============================================================================
# Analysis
# ==============================================================================
loginfo("Calculating")

# get differential expression status for each gene in each event
exprs[, DGE_Status := "None"]
exprs[(qval < 0.05) & (log2fc > 0), DGE_Status := "Up"]
exprs[(qval < 0.05) & (log2fc < 0), DGE_Status := "Down"]
exprs[, DGE_Status := factor(DGE_Status, levels = c("Down", "None", "Up"), ordered = TRUE)]

# merge this information with contacts
# first, merge the expression status with the original gene IDs
contacts_exprs <- merge(
    x = contacts,
    y = exprs[, .(event_ID, target_id, query_gene_dge_status = DGE_Status)],
    by.x = c("contact_event_ID", "query_gene_id"),
    by.y = c("event_ID", "target_id"),
    all.x = TRUE
)
# second, merge the expression status of the enhancer target genes
contacts_exprs <- merge(
    x = contacts_exprs,
    y = exprs[, .(event_ID, target_id, enh_target_gene_dge_status = DGE_Status)],
    by.x = c("contact_event_ID", "enh_target_gene_id"),
    by.y = c("event_ID", "target_id")
)

# state whether this contact is between the enhancer and
# 1. its original target gene
contacts_exprs[
    (query_gene_id == enh_target_gene_id),
    Contact_Type := "Self"
]
# 2. non-target genes that are over-expressed in the mutant
# 3. non-target genes that are under-expressed in the mutant
contacts_exprs[
    (query_gene_id != enh_target_gene_id),
    Contact_Type := paste("Other", enh_target_gene_dge_status)
]

# convert to a factor for easy plotting
contacts_exprs[, Contact_Type := factor(
    Contact_Type,
    levels = rev(c("Self", "Other Down", "Other None", "Other Up")),
    ordered = TRUE
)]

# calculate mean contacts between
contacts_exprs[, .N, by = c("query_gene_dge_status", "Contact_Type")]
agg_contacts <- contacts_exprs[,
    .(Mean_Obs_Exp_Freq = mean(Mean_Obs_Exp_Freq, na.rm = TRUE)),
    by = c("query_gene_dge_status", "Contact_Type", "sample_mutated_in_event", "query_gene_id", "query_gene_name")
]

# ==============================================================================
# Plots
# ==============================================================================
loginfo("Plotting figures")

gg <- (
    ggplot(
        data = contacts,
        mapping = aes(x = log2(Mean_Obs_Exp_Freq))
    )
    + geom_density()
    + geom_path(
        data = data.table(
            x = seq(-3, 3, 0.01),
            # y = dnorm(seq(-3, 3, 0.01), sd = 0.7270956)
            y = dnorm(seq(-3, 3, 0.01), sd = 0.6)
        ),
        mapping = aes(x = x, y = y)
    )
    + scale_x_continuous(
        # limits = c(0, 5)
    )
    + theme_minimal()
)
ggsave(
    file.path(PLOT_DIR, "dist.png"),
    gg,
    width = 12,
    height = 8,
    units = "cm"
)

gg_agg <- (
    ggplot(
        data = agg_contacts,
        mapping = aes(
            x = sample_mutated_in_event,
            y = log2(Mean_Obs_Exp_Freq)
        )
    )
    # + geom_point(
    #     position = position_jitter(height = 0)
    # )
    + geom_path(
        aes(group = query_gene_id)
    )
    + facet_grid(Contact_Type ~ query_gene_dge_status)
    + theme_minimal()
)
ggsave(
    file.path(PLOT_DIR, "agg-contacts.png"),
    gg_agg,
    width = 20,
    height = 20,
    units = "cm"
)
