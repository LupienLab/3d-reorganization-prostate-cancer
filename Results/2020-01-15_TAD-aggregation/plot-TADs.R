# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))

# ==============================================================================
# Functions
# ==============================================================================

# ==============================================================================
# Data
# ==============================================================================
# load metadata
metadata = fread(file.path("..", "..", "Data", "External", "LowC_Samples_Data_Available.tsv"))
SAMPLES = paste0("PCa", metadata[, get("Sample ID")])

# load aggregated boundary calls from each sample
boundaries = rbindlist(lapply(
    SAMPLES,
    function(s) {
        dt = fread(
            file.path("resolved-TADS", paste0(s, ".40000bp.aggregated-boundaries.tsv")),
            sep = "\t",
            header = TRUE,
            drop = "w"
        )
        dt[, SampleID := s]
        return(dt)
    }
))

# ==============================================================================
# Analysis
# ==============================================================================
boundary_counts = boundaries[, .N, by = "SampleID"]
boundary_counts_order = boundaries[, .N, by = c("SampleID", "Order")]

# ==============================================================================
# Plots
# ==============================================================================
# plot number of resolved boundaries
gg = (
    ggplot(data = boundary_counts)
    + geom_col(aes(x = SampleID, y = N, fill = SampleID))
    + scale_fill_viridis_d()
    + labs(x = NULL, y = "Number of unique boundaries")
    + guides(fill = FALSE)
    + theme_minimal()
    + theme(
        axis.text.x = element_text(angle = 90)
    )
)
ggsave(
    "Plots/tad-counts.png",
    height = 12,
    width = 20,
    units = "cm"
)

# plot number of resolved boundaries by order
gg = (
    ggplot(data = boundary_counts_order)
    + geom_col(aes(x = Order, y = N, fill = SampleID, group = SampleID), position = "dodge")
    + scale_fill_viridis_d()
    + xlab(c(0, 28))
    + labs(x = "Boundary Order", y = "Number of unique boundaries")
    + theme_minimal()
    + theme(
        axis.text.x = element_text(angle = 90)
    )
)
ggsave(
    "Plots/tad-counts-by-order.png",
    height = 12,
    width = 20,
    units = "cm"
)
