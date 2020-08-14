# ==============================================================================
# Meta
# ==============================================================================
# Plot BPScore
# --------------------------------------
# Description: Plot TAD similarities, calculated by BPscore
# Author: James Hawley


# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
suppressMessages(library("MASS"))
suppressMessages(library("pheatmap"))
source("savefig.R")

MAX_WINDOW <- 24
MIN_WINDOW <- 3
MAX_PERSISTENCE <- MAX_WINDOW - MIN_WINDOW + 1
PLOT_DIR <- "Plots"


# ==============================================================================
# Data
# ==============================================================================
cat("Loading Data\n")
# load metadata
metadata <- fread("config.tsv", sep = "\t", header = TRUE)
metadata <- metadata[Include == "Yes"]
SAMPLES <- metadata[, SampleID]
TUMOUR_SAMPLES <- metadata[Source == "Primary" & Type == "Malignant", SampleID]
LINE_SAMPLES <- metadata[Source == "Cell Line", SampleID]

# load BPscore calculations
bpscore <- fread("Statistics/tad-distances.tsv", sep = "\t", header = TRUE)
window_diffs <- fread("Statistics/tad-similarity-deltas.tsv", sep = "\t", header = TRUE)


# ==============================================================================
# Analysis
# ==============================================================================
# colour bpscore rows by primary-primary, line-line, or primary-line comparisons
cat("Counting similarity\n")
bpscore <- merge(
    bpscore,
    metadata[, .(SampleID, Source, Label, Type_Colour)],
    by.x = "s2",
    by.y = "SampleID",
    all.x = FALSE,
    all.y = TRUE
)
bpscore <- merge(
    bpscore,
    metadata[, .(SampleID, Source, Label, Type_Colour)],
    by.x = "s1",
    by.y = "SampleID",
    suffixes = c("_2", "_1"),
    all.x = FALSE,
    all.y = TRUE,
)
# if s1 and s2 come from the same source material, return Type_Colour, otherwise return green colour
bpscore[Source_1 == Source_2 & Source_1 == "Primary", Colour := "#1F77B4"]
bpscore[Source_1 == Source_2 & Source_1 == "Cell Line", Colour := "#FF7F0D"]
bpscore[Source_1 != Source_2, Colour := "#3FE686"]
bpscore[, `:=`(
    Cell_Type_1 = gsub(" .*", "", Label_1),
    Cell_Type_2 = gsub(" .*", "", Label_2)
)]
bpscore[, Cell_Type_Combo := paste(Cell_Type_1, Cell_Type_2, sep = "+")]

# calculate median similarity between different batches
non_prostate <- 1:6
prostate_primary <- 7:18
prostate_line <- 19:25

# cluster samples by BPscore distances
# recast data to wide form for plotting with pheatmap
mat <- as.matrix(
    dcast(bpscore[w == MAX_WINDOW], s1 ~ s2, value.var = "dist", fill = 0),
    rownames = "s1"
)

# compare primary prostate to non-prostate lines
cat("Non-prostate + Primary Prostate\n")
summary(as.vector(1 - mat[non_prostate, prostate_primary]))

# compare primary prostate to prostate lines
cat("Primary Prostate + Prostate Line\n")
summary(as.vector(1 - mat[prostate_primary, prostate_line]))

# compare non-prostate lines to prostate lines
cat("Non-prostate + Prostate Line\n")
summary(as.vector(1 - mat[non_prostate, prostate_line]))

# ==============================================================================
# Plots
# ==============================================================================
cat("Plotting\n")
# plot distances between TADs of different samples across window sizes
gg_bpscore_primary = (
    ggplot(data = bpscore[s1 != s2 & w <= MAX_WINDOW & !grepl("(4DN|SRR|BP)", s1) & !grepl("(4DN|SRR|BP)", s2)])
    + geom_point(
        # aes(x = w, y = 1 - dist, colour = Colour),
        aes(x = w, y = 1 - dist, colour = w),
        position = position_jitter(width = 0.2, height = 0)
    )
    + geom_boxplot(
        aes(x = w, y = 1 - dist, group = w),
        alpha = 0.5,
        outlier.shape = NA,
        width = 0.5
    )
    + labs(x = "Window size", y = "TAD similarity (1 - BPscore)")
    + scale_x_discrete(
        breaks = seq(MIN_WINDOW, MAX_WINDOW, by = 3),
        limits = seq(MIN_WINDOW, MAX_WINDOW, by = 3),
        labels = seq(MIN_WINDOW, MAX_WINDOW, by = 3)
    )
    + scale_colour_viridis_c()
    + guides(colour = FALSE)
    + theme_minimal()
)
savefig(gg_bpscore_primary, file.path(PLOT_DIR, "bp-score"))

# plot distances between TADs of different samples across window sizes
# (same as above, but including all samples, not just the primaries)
gg_bpscore_prostate_line = (
    ggplot(
        data = bpscore[s1 < s2 & w <= MAX_WINDOW & grepl("SRR", s1) & grepl("SRR", s2)],
        mapping = aes(x = w, y = 1 - dist, colour = Cell_Type_Combo)
    )
    + geom_point(
        position = position_jitter(width = 0.2, height = 0)
    )
    + geom_smooth(
        aes(group = Cell_Type_Combo, fill = Cell_Type_Combo),
        method = "loess",
        alpha = 0.5
    )
    + labs(x = "Window size", y = "TAD similarity (1 - BPscore)")
    + scale_x_discrete(
        breaks = seq(MIN_WINDOW, MAX_WINDOW, by = 3),
        limits = seq(MIN_WINDOW, MAX_WINDOW, by = 3),
        labels = seq(MIN_WINDOW, MAX_WINDOW, by = 3)
    )
    + guides(colour = guide_legend(title = "Comparison"), fill = FALSE)
    + theme_minimal()
)
savefig(gg_bpscore_prostate_line, file.path(PLOT_DIR, "bp-score.prostate-lines"))

# plot the change in similarities over window sizes
gg_sim_window = (
    ggplot(data = window_diffs[w <= MAX_WINDOW])
    + geom_path(aes(x = w, y = 1 - diff, group = SampleID, colour = SampleID))
    + labs(x = "Window size", y = expression("1 - " * delta["w, w-1"]))
    + scale_x_discrete(
        breaks = seq(MIN_WINDOW, MAX_WINDOW, by = 3),
        limits = seq(MIN_WINDOW, MAX_WINDOW, by = 3),
        labels = seq(MIN_WINDOW, MAX_WINDOW, by = 3)
    )
    + scale_colour_manual(
        limits = metadata[, SampleID],
        labels = metadata[, Label],
        values = metadata[, Sample_Colour],
        name = "Sample"
    )
    + theme_minimal()
)
savefig(gg_sim_window, file.path(PLOT_DIR, "bp-score.window-similarity"))


# plot the change in similarities over window sizes
gg_sim_window_delta = (
    ggplot(data = window_diffs[w >= 5 & w <= MAX_WINDOW])
    + geom_path(aes(x = w, y = abs_delta, group = SampleID, colour = SampleID))
    + labs(x = "Window size", y = expression("|" * delta["w, w-1"] - delta["w-1, w-2"] * "|"))
    + scale_x_discrete(
        breaks = seq(MIN_WINDOW, MAX_WINDOW, by = 3),
        limits = seq(MIN_WINDOW, MAX_WINDOW, by = 3),
        labels = seq(MIN_WINDOW, MAX_WINDOW, by = 3)
    )
    + scale_colour_manual(
        limits = metadata[, SampleID],
        labels = metadata[, Label],
        values = metadata[, Sample_Colour],
        name = "Patient"
    )
    + theme_minimal()
)
savefig(gg_sim_window_delta, file.path(PLOT_DIR, "bp-score.window-similarity.delta"))

# heatmap of similarity at MAX_WINDOW
annot_col <- as.data.frame(metadata[, .SD, .SDcols = c("Type", "Tissue", "Source")])
rownames(annot_col) <- metadata[, SampleID]

for (ext in c("png", "pdf")) {
    pheatmap(
        mat = 1 - mat,
        filename = paste0(file.path(PLOT_DIR, "bp-score.cluster."), ext),
        clustering_rows = TRUE,
        clustering_cols = TRUE,
        clustering_distance_rows = as.dist(mat),
        clustering_distance_cols = as.dist(mat),
        clustering_method = "complete",
        labels_row = metadata[order(SampleID), Label],
        labels_col = metadata[order(SampleID), Label],
        legend = TRUE,
        show_rownames = FALSE,
        annotation_col = annot_col
    )
}
