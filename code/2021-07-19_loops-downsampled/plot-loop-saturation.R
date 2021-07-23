# ==============================================================================
# Meta
# ==============================================================================
# loop-saturation
# --------------------------------------
# Description: Plot the loop call saturation estimates
# Author: James Hawley

# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("logging"))

loginfo("Loading packages")
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
suppressMessages(library("matrixStats"))
source(file.path("..", "src", "savefig.R"))
source(file.path("..", "src", "theme.R"))

set.seed(42)

RES_DIR <- file.path("..", "..", "results", "2021-07-19_loops-downsampled")
LOOP_DIR <- file.path(RES_DIR, "Loops")
SAT_DIR <- file.path(RES_DIR, "Saturation")
PLOT_DIR <- file.path(RES_DIR, "Plots")


# ==============================================================================
# Data
# ==============================================================================
loginfo("Loading data")

metadata <- fread("config.tsv", sep = "\t", header = TRUE)
metadata <- metadata[Include == "Yes", .SD]

saturation_models <- readRDS(
    file.path(SAT_DIR, "loop-saturation.models.rds")
)

all_loop_counts <- fread(
    file.path(SAT_DIR, "loop-saturation.iterations.tsv"),
    sep = "\t",
    header = TRUE
)

depth_ests <- fread(
    file.path(SAT_DIR, "loop-saturation.depth-estimates.tsv"),
    sep = "\t",
    header = TRUE
)

saturation_ests <- fread(
    file.path(SAT_DIR, "loop-saturation.saturation-estimates.tsv"),
    sep = "\t",
    header = TRUE
)

obs_saturation_ests <- fread(
    file.path(SAT_DIR, "loop-saturation.observed-saturation.tsv"),
    sep = "\t",
    header = TRUE
)

asymptotes <- rbindlist(lapply(
    names(saturation_models),
    function(s) {
        data.table(
            SampleID = s,
            Label = metadata[SampleID == s, Label],
            Seq_Depth = 1e9,
            N = depth_ests[(SampleID == s) & (Seq_Depth == 1e9), Est_N_Loops]
        )
    }
))

# ==============================================================================
# Plots
# ==============================================================================
loginfo("Plotting data")

for (s in names(saturation_models)) {
    gg <- (
        ggplot()
        + geom_point(
            data = all_loop_counts[SampleID == s],
            mapping = aes(x = Seq_Depth, y = N, colour = Seq_Depth),
            position = position_jitter(width = 0.2, height = 0),
            alpha = 0.1
        )
        + geom_boxplot(
            data = all_loop_counts[SampleID == s],
            mapping = aes(x = Seq_Depth, y = N, group = Seq_Depth, fill = Seq_Depth),
            outlier.shape = NA,
            alpha = 0.5
        )
        # asymptotic curve fit
        + geom_path(
            data = depth_ests[SampleID == s],
            aes(x = Seq_Depth, y = Est_N_Loops),
            linetype = "dashed"
        )
        # label for sample
        + geom_text(
            aes(
                x = saturation_models[[ s ]]$estimates$saturation[Frac_Saturation >= 0.95, mean(Seq_Depth)],
                y = 0.9 * saturation_models[[ s ]]$model$Asym,
                label = "Model Fit"
            ),
            vjust = 0,
            hjust = 0.5
        )
        # add asymptotic value
        + geom_hline(
            aes(yintercept = c(
                saturation_models[[ s ]]$model$Asym,
                saturation_models[[ s ]]$interactions$identified
            )),
            linetype = "dashed",
            colour = "#1e90ff"

        )
        # labels for horizontal lines
        + geom_text(
            aes(
                x = c(0, 0),
                y = c(
                    saturation_models[[ s ]]$interactions$identified,
                    saturation_models[[ s ]]$model$Asym
                ),
                label = c(
                    paste0(
                        "Interactions Detected: ",
                        saturation_models[[ s ]]$interactions$identified,
                        " (", 100 * round(saturation_models[[ s ]]$interactions$saturation, 3),
                        "%)"
                    ),
                    paste(
                        "Estimated Number of Interactions:",
                        round(saturation_models[[ s ]]$model$Asym)
                    )
                )
            ),
            vjust = -1,
            hjust = 0,
            colour = "#1e90ff"
        )
        # add bars for number of samples required to reach saturation
        + geom_vline(
            aes(xintercept = saturation_models[[ s ]]$estimates$saturation[, round(Seq_Depth)]),
            linetype = "dashed",
            colour = "#66cdaa"
        )
        # labels for vertical lines
        + geom_text(
            data = saturation_models[[ s ]]$estimates$saturation,
            aes(
                x = round(Seq_Depth),
                y = 0,
                label = paste0(
                    100 * Frac_Saturation,
                    "%: ",
                    round(Seq_Depth)
                )
            ),
            vjust = -0.5,
            hjust = 0,
            angle = 90,
            colour = "#66cdaa"
        )
        + scale_fill_viridis_c()
        + scale_colour_viridis_c()
        + scale_x_continuous(
            name = "Sequencing Depth (M)",
            limits = c(
                0,
                # estimated number of samples required to reach 99% saturation of loop calls
                saturation_models[[ s ]]$estimates$saturation[
                    Frac_Saturation == 0.99,
                    round(Seq_Depth) + 1
                ]
            ),
            breaks = seq(500000000, 2e9, 500000000),
            labels = seq(500, 2e3, 500)
        )
        + scale_y_continuous(
            name = "Number of interactions",
            limits = c(0, 1.05 * saturation_models[[ s ]]$model$Asym)
        )
        + guides(fill = FALSE, colour = FALSE)
        + theme_minimal()
        + jrh_theme()
    )
    savefig(gg, file.path(PLOT_DIR, paste0("loop-saturation.", s)))
}

gg <- (
    ggplot()
    # asymptotic curve fit
    + geom_path(
        data = depth_ests,
        mapping = aes(x = Seq_Depth, y = Est_N_Loops, colour = SampleID),
        linetype = "dashed"
    )
    # boxplots of downsampling iterations
    + geom_boxplot(
        data = all_loop_counts[!is.na(Replicate)],
        mapping = aes(
            x = Seq_Depth,
            y = N,
            group = Seq_Depth
        ),
        colour = "#000000",
        alpha = 0.1
    )
    # points of downsampling iterations
    + geom_point(
        data = all_loop_counts[!is.na(Replicate)],
        mapping = aes(
            x = Seq_Depth,
            y = N,
            fill = SampleID
        ),
        colour = "#000000",
        shape = 21,
        size = 4,
        alpha = 0.05
    )
    # points of observed sequencing depths and interaction counts
    + geom_point(
        data = all_loop_counts[is.na(Replicate)],
        mapping = aes(
            x = Seq_Depth,
            y = N,
            fill = SampleID
        ),
        colour = "#000000",
        shape = 21,
        size = 4
    )
    + scale_x_continuous(
        name = "Sequencing Depth"
    )
    + scale_y_continuous(
        name = "Number of interactions",
        sec.axis = dup_axis(
            breaks = pmin(asymptotes$N, 20000),
            labels = asymptotes$Label,
            name = NULL
        )
    )
    + guides(colour = "none", fill = "none")
    + theme_minimal()
    + jrh_theme()
)
savefig(gg, file.path(PLOT_DIR, "loop-saturation.all"))
