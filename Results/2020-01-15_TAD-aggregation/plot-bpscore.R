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
    metadata[, .(SampleID, Source, Tissue, Label, Type)],
    by.x = "s1",
    by.y = "SampleID",
    all.x = FALSE,
    all.y = TRUE
)
bpscore <- merge(
    bpscore,
    metadata[, .(SampleID, Source, Tissue, Label, Type)],
    by.x = "s2",
    by.y = "SampleID",
    suffixes = c("_1", "_2"),
    all.x = FALSE,
    all.y = TRUE,
)

# sample comparisons for various plots
bpscore[(Tissue_1 == "Prostate") & (Tissue_2 == "Prostate"), Prostate_Type_Combo := paste(Type_1, Type_2, sep = "+")]
bpscore[
    (Tissue_1 == "Prostate") & (Tissue_2 == "Prostate") & (Source_1 == "Cell Line") & (Source_2 == "Cell Line"),
    Cell_Type_Combo := paste(
        gsub(" .*", "", Label_1),
        gsub(" .*", "", Label_2),
        sep = "+"
    )
]

# hypothesis testing for differences between TADs
# (similarity of Benign-vs-Benign comparisons to Benign-vs-Tumour)
# (and likewise, Tumour-vs-Tumour to Benign-vs-Tumour)
bpscore_primary_tumour_vs_benign <- bpscore[
    s1 < s2
    & w == MAX_WINDOW
    & (Source_1 == "Primary" & Source_2 == "Primary")
]

htests <- list(
    "primary_benign-only_benign-vs-tumour" = wilcox.test(
        x = bpscore_primary_tumour_vs_benign[Prostate_Type_Combo == "Benign+Benign", 1 - dist],
        y = bpscore_primary_tumour_vs_benign[Prostate_Type_Combo == "Benign+Malignant", 1 - dist],
        paired = FALSE
    ),
    "primary_tumour-only_benign-vs-tumour" = wilcox.test(
        x = bpscore_primary_tumour_vs_benign[Prostate_Type_Combo == "Malignant+Malignant", 1 - dist],
        y = bpscore_primary_tumour_vs_benign[Prostate_Type_Combo == "Benign+Malignant", 1 - dist],
        paired = FALSE
    )
)

# cluster samples by BPscore distances
# recast data to wide form for plotting with pheatmap
mat <- as.matrix(
    dcast(bpscore[w == MAX_WINDOW], s1 ~ s2, value.var = "dist", fill = 0),
    rownames = "s1"
)

# compare different groups of samples
group_idx <- list(
    "Primary prostate tumour" = 12:23,
    "Primary prostate benign" = 7:11,
    "Cell line prostate tumour" = c(24:25, 28:30),
    "Cell line prostate benign" = 26:27,
    "Cell line non-prostate" = 1:6
)
median_mat <- matrix(ncol = length(group_idx), nrow = length(group_idx))
for (i in 1:length(group_idx)) {
    for (j in i:length(group_idx)) {
        median_mat[i, j] <- median(as.vector(1 - mat[group_idx[[i]], group_idx[[j]]]))
    }
}
rownames(median_mat) <- names(group_idx)
colnames(median_mat) <- rownames(median_mat)

cat("Sample type comparisons\n")
print(median_mat)
write.matrix(
    median_mat,
    file.path("Statistics", "tad-similarity.sample-types.tsv"),
    sep = "\t"
)

# ==============================================================================
# Plots
# ==============================================================================
cat("Plotting\n")
# plot distances between TADs of different primary tumour samples across window sizes
gg_bpscore_primary_tumour = (
    ggplot(data = bpscore[s1 < s2 & w <= MAX_WINDOW & grepl("PCa", s1) & grepl("PCa", s2)])
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
savefig(gg_bpscore_primary_tumour, file.path(PLOT_DIR, "bp-score.primary.tumour"))

# plot distances between TADs of different primary tumour samples across window sizes
gg_bpscore_primary_benign = (
    ggplot(data = bpscore[s1 < s2 & w <= MAX_WINDOW & grepl("Benign-Prostate", s1) & grepl("Benign-Prostate", s2)])
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
savefig(gg_bpscore_primary_benign, file.path(PLOT_DIR, "bp-score.primary.benign"))

# primary samples: benign vs tumour
gg_bpscore_primary_tumour_vs_benign <- (
    ggplot()
    + geom_smooth(
        data = bpscore[
            s1 < s2
            & w <= MAX_WINDOW
            & (Source_1 == "Primary" & Source_2 == "Primary")
            & (Source_1 == "Primary" & Source_2 == "Primary")
        ],
        mapping = aes(
            x = w,
            y = 1 - dist,
            colour = Prostate_Type_Combo,
            fill = Prostate_Type_Combo,
            group = Prostate_Type_Combo
        ),
        method = "loess",
        alpha = 0.5
    )
    + labs(x = "Window size", y = "TAD similarity (1 - BPscore)")
    + scale_x_discrete(
        breaks = seq(MIN_WINDOW, MAX_WINDOW, by = 3),
        limits = seq(MIN_WINDOW, MAX_WINDOW, by = 3),
        labels = seq(MIN_WINDOW, MAX_WINDOW, by = 3)
    )
    + scale_colour_discrete(
        name = "Comparison",
        breaks = c("Benign+Benign", "Benign+Malignant", "Malignant+Malignant"),
        labels = c("Benign vs Benign", "Benign vs Tumour", "Tumour vs Tumour")
    )
    + scale_fill_discrete(
        name = "Comparison",
        breaks = c("Benign+Benign", "Benign+Malignant", "Malignant+Malignant"),
        labels = c("Benign vs Benign", "Benign vs Tumour", "Tumour vs Tumour")
    )
    + guides(fill = FALSE)
    + theme_minimal()
    + theme(legend.position = "bottom")
)
savefig(gg_bpscore_primary_tumour_vs_benign, file.path(PLOT_DIR, "bp-score.primary.tumour-vs-benign"))

# same as above, but using the sample standard deviation for the ribbon
# instead of the 95% confidence interval around the regressed mean
gg_bpscore_primary_tumour_vs_benign <- (
    ggplot()
    + geom_ribbon(
        data = bpscore[
            s1 < s2
            & w <= MAX_WINDOW
            & (Source_1 == "Primary" & Source_2 == "Primary")
            & (Source_1 == "Primary" & Source_2 == "Primary"),
            # using `get` because `dist` is a builtin function name
            .(Mean = mean(1 - dist), SD = sd(1 - dist)),
            by = c("w", "Prostate_Type_Combo")
        ],
        mapping = aes(
            x = w,
            ymin = Mean - SD,
            ymax = Mean + SD,
            fill = Prostate_Type_Combo
        ),
        alpha = 0.2
    )
    + geom_smooth(
        data = bpscore[
            s1 < s2
            & w <= MAX_WINDOW
            & (Source_1 == "Primary" & Source_2 == "Primary")
            & (Source_1 == "Primary" & Source_2 == "Primary")
        ],
        mapping = aes(
            x = w,
            y = 1 - dist,
            colour = Prostate_Type_Combo,
            fill = Prostate_Type_Combo,
            group = Prostate_Type_Combo
        ),
        method = "loess",
        se = FALSE,
        alpha = 0.5
    )
    + labs(x = "Window size", y = "TAD similarity (1 - BPscore)")
    + scale_x_discrete(
        breaks = seq(MIN_WINDOW, MAX_WINDOW, by = 3),
        limits = seq(MIN_WINDOW, MAX_WINDOW, by = 3),
        labels = seq(MIN_WINDOW, MAX_WINDOW, by = 3)
    )
    + scale_colour_discrete(
        name = "Comparison",
        breaks = c("Benign+Benign", "Benign+Malignant", "Malignant+Malignant"),
        labels = c("Benign vs Benign", "Benign vs Tumour", "Tumour vs Tumour")
    )
    + scale_fill_discrete(
        name = "Comparison",
        breaks = c("Benign+Benign", "Benign+Malignant", "Malignant+Malignant"),
        labels = c("Benign vs Benign", "Benign vs Tumour", "Tumour vs Tumour")
    )
    + guides(fill = FALSE)
    + theme_minimal()
    + theme(legend.position = "bottom")
)
savefig(gg_bpscore_primary_tumour_vs_benign, file.path(PLOT_DIR, "bp-score.primary.tumour-vs-benign.stdev"))

# cell lines: benign vs tumour
gg_bpscore_line_tumour_vs_benign = (
    ggplot(
        data = bpscore[
            s1 < s2
            & w <= MAX_WINDOW
            & (Source_1 == "Cell Line" & Source_2 == "Cell Line")
            & (Tissue_1 == "Prostate" & Tissue_2 == "Prostate")
        ],
        mapping = aes(
            x = w,
            y = 1 - dist,
            colour = Cell_Type_Combo,
            fill = Cell_Type_Combo,
            group = Cell_Type_Combo
        )
    )
    + geom_smooth(
        method = "loess",
        alpha = 0.5
    )
    + labs(x = "Window size", y = "TAD similarity (1 - BPscore)")
    + scale_x_discrete(
        breaks = seq(MIN_WINDOW, MAX_WINDOW, by = 3),
        limits = seq(MIN_WINDOW, MAX_WINDOW, by = 3),
        labels = seq(MIN_WINDOW, MAX_WINDOW, by = 3)
    )
    + scale_colour_discrete(
        name = "Comparison"
    )
    + guides(fill = FALSE)
    + theme_minimal()
    + theme(legend.position = "bottom")
)
savefig(gg_bpscore_line_tumour_vs_benign, file.path(PLOT_DIR, "bp-score.line.tumour-vs-benign"))


# plot the change in similarities over window sizes
gg_sim_window = (
    ggplot(data = window_diffs[w <= MAX_WINDOW])
    + geom_path(aes(x = w, y = 1 - BPscore, group = SampleID, colour = SampleID))
    + labs(x = "Window size", y = "Similarity between consecutive window sizes\n(1 - BPscore)")
    + scale_x_discrete(
        breaks = seq(MIN_WINDOW, MAX_WINDOW, by = 3),
        limits = seq(MIN_WINDOW, MAX_WINDOW, by = 3),
        labels = seq(MIN_WINDOW, MAX_WINDOW, by = 3)
    )
    + scale_colour_discrete(
        limits = metadata[, SampleID],
        labels = metadata[, Label],
        name = "Sample"
    )
    + theme_minimal()
)
savefig(gg_sim_window, file.path(PLOT_DIR, "bp-score.window-similarity"))


# plot the change in similarities over window sizes
gg_sim_window_delta = (
    ggplot(data = window_diffs[w >= 5 & w <= MAX_WINDOW])
    + geom_path(aes(x = w, y = window_step_abs_delta, group = SampleID, colour = SampleID))
    + labs(x = "Window size", y = "Gain of step-wise similarity")
    + scale_x_discrete(
        breaks = seq(MIN_WINDOW, MAX_WINDOW, by = 3),
        limits = seq(MIN_WINDOW, MAX_WINDOW, by = 3),
        labels = seq(MIN_WINDOW, MAX_WINDOW, by = 3)
    )
    + scale_colour_discrete(
        limits = metadata[, SampleID],
        labels = metadata[, Label],
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
        clustering_method = "ward.D2",
        treeheight_col = 0,
        labels_row = metadata[order(SampleID), Label],
        labels_col = metadata[order(SampleID), Label],
        legend = TRUE,
        show_rownames = FALSE,
        annotation_col = annot_col
    )
}
