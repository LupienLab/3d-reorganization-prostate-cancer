# ==============================================================================
# Meta
# ==============================================================================
# plot-sv-changes
# --------------------------------------
# Description: Summary statistics of the SVs and how they change gene expression
# Author: James Hawley


# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("logging"))

loginfo("Loading packages")
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
source(file.path("..", "src", "savefig.R"))

RES_DIR <- file.path(
    "..", "..", "results", "2020-06-18_sv-disruption-expression"
)
PLOT_DIR <- file.path(RES_DIR, "Plots")

# ==============================================================================
# Data
# ==============================================================================
loginfo("Loading data")

sv_disruption <- fread("summary-sv-disruption-annotation.tsv")


# ==============================================================================
# Plots
# ==============================================================================
loginfo("Plotting figures")

gg <- (
    ggplot(data = sv_disruption)
    + geom_bar(
        aes(x = status, fill = status),
        colour = "black"
    )
    + guides(fill = FALSE)
    + scale_x_discrete(
        name = "Expression Changes from SV",
        breaks = c("down only", "up and down", "up only"),
        labels = c("Under-Expression Only", "Over- and Under-Expression", "Over-Expression Only")
    )
    + scale_fill_manual(
        name = "Count",
        breaks = c("down only", "up and down", "up only"),
        labels = c("Under-Expression Only", "Over- and Under-Expression", "Over-Expression Only"),
        values = c("#4139E1", "#F666E9", "#FF6347")
    )
    + theme_minimal()
)
savefig(gg, "Plots/event-changes.summary", width = 7, height = 10)

gg <- (
    ggplot(data = sv_disruption[status == "up and down"])
    + geom_bar(
        aes(x = opposite_end_of_break, fill = direct_contact)
    )
    + scale_x_discrete(
        name = "Differential Expression On Opposite Ends of SV",
        breaks = c("yes", "no"),
        labels = c("Yes", "No")
    )
    + scale_y_continuous(
        name = "Number of SVs",
        limit = c(0, 15),
        breaks = c(0, 5, 10, 15)
    )
    + scale_fill_manual(
        name = "Contact Between Genes",
        breaks = c("yes", "no"),
        labels = c("Direct", "Indirect"),
        values = c("#4139E1", "#AEA28E")
    )
    + theme_minimal()
    + theme(
        legend.position = "bottom"
    )
)
savefig(gg, "Plots/event-changes.contact", width = 7, height = 10)
