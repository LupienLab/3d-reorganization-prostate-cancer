# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
suppressMessages(library("grid"))
suppressMessages(library("gridExtra"))
suppressMessages(library("gtable"))

# ==============================================================================
# Data
# ==============================================================================
# read sample metadata
cpcgene_metadata = fread(
    "../../Data/External/CPC-GENE/CPC-GENE_Fraser-2017_Nature_Clinical-Data.tsv",
    sep = "\t",
    header = TRUE
)

# remove `-F1` from `SampleID` column to match the values in `lowc_metadata`
cpcgene_metadata[, SampleID := gsub("-F1$", "", SampleID)]

# read Low-C sample IDs and metadata
lowc_metadata = fread(
    "../../Data/External/LowC_Samples_Data_Available.tsv",
    sep = "\t",
    header = TRUE
)
lowc_metadata <- metadata[Include == "Yes"]
# change column names for easier access
colnames(lowc_metadata)[1:2] = c("PCaID", "SampleID")
lowc_metadata[, PCaID := paste0("PCa", PCaID)]

# merge two metadata files together
metadata = merge(
    x = cpcgene_metadata,
    y = lowc_metadata,
    by = "SampleID",
    all.x = FALSE,
    all.y = TRUE
)
colnames(metadata)[18] = "ERG"

# read in breakpoints called by breakfinder
breakpoints = fread(
    "Breakpoints/Default/breakpoints.tsv",
    sep = "\t",
    header = TRUE
)

# ==============================================================================
# Analysis
# ==============================================================================
# count total number of breakpoints
breakpoints_counted = breakpoints[, .N, by = "Sample"][order(-N)]
colnames(breakpoints_counted) = c("PCaID", "Breakpoints")

# metadata grid for plotting
metadata_wide = melt(
    metadata,
    id.vars = "PCaID",
    measure.vars = c("Gleason Score", "BCR", "ERG"),
    variable.name = "Feature",
    value.name = "Value"
)
# order by sample in the same order as breakpoints_counted
metadata_wide[, PCaID := factor(PCaID, levels = breakpoints_counted$PCaID, ordered = TRUE)]
metadata_wide = metadata_wide[order(PCaID)]

# ==============================================================================
# Plots
# ==============================================================================
# plot according to number of SVs
gg_bars = (
    ggplot(data = breakpoints_counted)
    + geom_col(aes(x = factor(PCaID, ordered = TRUE, levels = PCaID), y = Breakpoints, fill = ))
    + labs(x = NULL, y = "SVs Detected")
    + theme_minimal()
    + theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank()
    )
)

# heatmap with sample metadata
gg_meta = (
    ggplot(data = metadata_wide)
    + geom_tile(aes(x = PCaID, y = Feature, fill = Value))
    + labs(x = "Sample", y = "Feature")
    + guides(fill = FALSE)
    + scale_fill_manual(
        limits = c("Negative", "Positive", "Yes", "No", "3+3", "3+4", "4+3"),
        values = c("#548235", "#2E75B6", "#000000", "#FFFFFF", "#EDB58F", "#ED985E", "#ED7D31")
        # name = "Colour",
        # labels = "",
        # breaks = 
    )
    + theme_minimal()
    + theme(
        axis.text.x = element_text(angle = 90),
        panel.grid.major = element_blank()
        # axis.text.y = element_blank()
    )
)

# combine plots
bars = ggplotGrob(gg_bars)
features = ggplotGrob(gg_meta)
g = rbind(bars, features, size = "first")
g$widths = unit.pmax(bars$widths, features$widths)

png("Plots/counted-breakpoints.png", width = 20, height = 12, units = "cm", res = 300)
grid.draw(g)
dev.off()
dev.off()
