# ==============================================================================
# Meta
# ==============================================================================
# CTCF enrichment
# --------------------------------------
# Description: Calculate and plot CTCF binding site enrichment around TAD boundaries
# Author: James Hawley


# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("argparse"))
source("savefig.R")

if (!interactive()) {
    PARSER <- argparse::ArgumentParser(
        description = "Calculate and plot CTCF binding site enrichment around TAD boundaries"
    )
    PARSER$add_argument(
        "cell",
        type = "character",
        help = "Cell line with CTCF binding sites to load",
        choices = c("22Rv1", "C4-2B", "LNCaP", "VCaP")
    )
    ARGS <- PARSER$parse_args()
} else {
    ARGS <- list(
        "cell" = "22Rv1"
    )
}

MAX_WINDOW <- 24
MIN_WINDOW <- 3
MAX_PERSISTENCE <- MAX_WINDOW - MIN_WINDOW + 1
PLOT_DIR <- "Plots"


# ==============================================================================
# Data
# ==============================================================================
cat("Loading data\n")
# load metadata
metadata <- fread("config.tsv", sep = "\t", header = TRUE)
metadata <- metadata[Include == "Yes"]
SAMPLES <- metadata[, SampleID]
TUMOUR_SAMPLES <- metadata[Source == "Primary" & Type == "Malignant", SampleID]
LINE_SAMPLES <- metadata[Source == "Cell Line", SampleID]

# load CTCF distances
ctcf_pairs = fread(
    file.path("CTCF", paste0("TAD-boundary.", ARGS$cell, "-CTCF-peaks.distances.tsv")),
    sep = "\t",
    header = TRUE
)
ctcf_pairs <- merge(ctcf_pairs, metadata[, .(SampleID, Source)])


# ==============================================================================
# Analysis
# ==============================================================================
cat("Calculating enrichment\n")
# calculate frequency at boundary (bin distance from boundary)
ctcf_fc <- ctcf_pairs[Bin_Mid == 0, .(Peak = Freq), keyby = c("SampleID", "Source")]

# take mean of enrichment >= 100 kbp away from nearest boundary
ctcf_fc$Background <- ctcf_pairs[abs(Bin_Mid) > 100000, mean(Freq), keyby = "SampleID"]$V1

# calculate fold change between distal background and proximal peak
ctcf_fc[, Fold := Peak / Background]
fwrite(
    ctcf_fc,
    file.path("Statistics", paste0("tad-boundary-enrichment.", ARGS$cell, "-CTCF.tsv")),
    sep = "\t",
    col.names = TRUE
)


# ==============================================================================
# Plots
# ==============================================================================
cat("Plotting\n")
# CTCF binding site proximity to boundaries
gg_bounds_ctcf <- (
    ggplot(data = ctcf_pairs[SampleID %in% c(TUMOUR_SAMPLES, LINE_SAMPLES)])
    + geom_path(aes(x = Bin_Mid / 1e3, y = Freq, colour = Source, group = SampleID))
    + labs(x = "Distance from TAD boundary (kbp)", y = paste("Average #", ARGS$cell, "CTCF Peaks / 5 kbp"))
    + scale_x_continuous(
        breaks = seq(-150, 150, 50),
        labels = seq(-150, 150, 50),
    )
    + scale_y_continuous(
        limits = c(0, 0.03)
    )
    + scale_colour_manual(
        limits = metadata[SampleID %in% c(TUMOUR_SAMPLES, LINE_SAMPLES), unique(Source)],
        labels = metadata[SampleID %in% c(TUMOUR_SAMPLES, LINE_SAMPLES), unique(Source)],
        values = c("#1F77B4", "#FF7F0D"),
        name = "Source"
    )
    + coord_cartesian(xlim = c(-150, 150))
    + theme_minimal()
)
savefig(gg_bounds_ctcf, file.path(PLOT_DIR, paste0("boundary-counts.ctcf-proximity.", ARGS$cell)))

gg_bounds_ctcf_fc <- (
    ggplot(data = ctcf_fc[SampleID %in% TUMOUR_SAMPLES | grepl("SRR", SampleID)])
    + geom_col(aes(x = SampleID, y = Fold, fill = SampleID))
    + labs(x = NULL, y = "Fold change (peak vs background)")
    + scale_x_discrete(
        breaks = metadata[, SampleID],
        labels = metadata[, Label]
    )
    + scale_fill_manual(
        limits = metadata[, SampleID],
        labels = metadata[, Label],
        values = metadata[, Type_Colour],
        name = "Patient"
    )
    + guides(fill = FALSE)
    + theme_minimal()
    + theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    )
)
savefig(gg_bounds_ctcf_fc, file.path(PLOT_DIR, paste0("boundary-counts.ctcf-proximity.", ARGS$cell, ".fold")))
