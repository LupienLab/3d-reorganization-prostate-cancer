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
colnames(metadata) = gsub(" ", "_", colnames(metadata))
metadata[, Sample_ID := paste0("PCa", Sample_ID)]

# load the aggregated counts for each sample
contact_dists = rbindlist(lapply(
    metadata[, Sample_ID],
    function(id) {
        dt = fread(file.path("Contacts", paste0(id, ".tsv")))
        colnames(dt) = c("Distance", "Count")
        dt[, SampleID := id]
        return(dt)
    }
))

# ==============================================================================
# Analysis
# ==============================================================================
mean_dists = contact_dists[, mean(Count), by = "Distance"]
colnames(mean_dists) = c("Distance", "Count")

# ==============================================================================
# Plots
# ==============================================================================
gg = (
    ggplot()
    # + geom_path(
    #     data = contact_dists,
    #     mapping = aes(x = Distance, y = Count, group = SampleID, colour = SampleID)
    # )
    + geom_smooth(
        data = mean_dists,
        mapping = aes(x = Distance, y = Count)
    )
    + scale_x_log10()
    + scale_y_log10()
    + labs(x = "Distance (bp)", y = expression("Count")
)
ggsave(
    file.path("Plots", "sample-contact-distance-distribution.png"),
    height = 12,
    width = 20,
    units = "cm"
)
