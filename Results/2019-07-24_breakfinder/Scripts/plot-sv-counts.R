# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))

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

# read in breakpoints called by breakfinder
breakpoints = fread(
    "Breakpoints/Default/breakpoints.tsv",
    sep = "\t",
    header = TRUE
)

# ==============================================================================
# Analysis
# ==============================================================================
breakpoints_counted = breakpoints[, .N, by = "Sample"][order(-N)]

# ==============================================================================
# Plots
# ==============================================================================
gg = (
    ggplot(data = breakpoints_counted)
    + geom_col(aes(x = Sample, y = N))
    + theme_minimal()
)
ggsave(
    "Plots/counted-breakpoints.png",
    width = 20,
    height = 12,
    units = "cm"
)
