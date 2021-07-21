# ==============================================================================
# Meta
# ==============================================================================
# loop-saturation
# --------------------------------------
# Description: Calculate and plot the loop call saturation estimates
# Author: James Hawley
# Steps:
#   1. On a cohort of samples, the sample order is randomised and the number of
#      peaks/regions identified in the first sample is plotted at x=1
#   2. The peaks identified in the second sample are compared to those found in
#      the first and only those which are new/unique are kept.
#      The sum of these new Sample2-only peaks and the Sample1-peaks is then
#      plotted for x=2
#   3. For every following sample, the number of NEW peaks identified in that
#      sample is added to the cumulative sum of unique peaks identified in all
#      preceding samples and plotted at x=sample number
#   4. Steps 1-3 are repeated n times. n should be defined by the user.
#   5. The final plot displays the cumulative peak number at each sample number
#      averaged over the n iterations, as well as standard error of the mean
#      error bars.
#   6. A self-starting non-linear regression model is fitted to the data and is
#      used to predict at what sample number saturation will be reached and what
#      level of saturation has been reached at the current sample size.


# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("logging"))

loginfo("Loading packages")
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
suppressMessages(library("matrixStats"))
source(file.path("..", "src", "savefig.R"))

set.seed(42)

RES_DIR <- file.path("..", "..", "results", "2021-07-19_loops-downsampled")
LOOP_DIR <- file.path(RES_DIR, "Loops")
SAT_DIR <- file.path(RES_DIR, "Saturation")
PLOT_DIR <- file.path(RES_DIR, "Plots")

DEPTHS <- as.integer(50000000 * (1:6))

# ==============================================================================
# Functions
# ==============================================================================
#' Fit exponential model to the loop counts to estimate the total number of loops
#' present at a given sequencing depth
#'
#' @param loop_counts Loop counts
#' @param max_depth maximum sequencing depth to estimate the model out to
#' @return structured list
estimate_saturation <- function(loop_counts, max_depth = 1e9) {
    # 1. Fit the exponential
    # ---------------------------------
    # use self starting model for asymptotic data SSasymp
    # it follows the formula: y = Asym + (R0 - Asym) * exp(-exp(lrc) * x) where
    # x     a numeric vector of values at which to evaluate the model
    # Asym  a numeric parameter for the upper asymptotic value of the model.
    # R0    a numeric parameter for the response when input is zero
    # lrc	a numeric parameter for the natural log of the rate constant
    nlsfitSS <- nls(
        N ~ SSasymp(Seq_Depth, Asym, R0, lrc),
        data = loop_counts
    )

    # predict the future values with increasing numbers of samples
    pred <- predict(
        nlsfitSS,
        list(N = seq(50000000, max_depth, 50000000))
    )

    # check that the model is a good fit for the data
    #   Residual sum of squares
    RSS <- sum(residuals(nlsfitSS)^2)
    #   Total sum of squares
    TSS <- loop_counts[, sum((N - mean(N))^2)]
    #   R-squared measure (this should be close to 1)
    r2 <- 1 - (RSS/TSS)

    # extract key parameters from the model
    Asym = summary(nlsfitSS)$coefficients[1]
    R0 = summary(nlsfitSS)$coefficients[2]
    lrc = summary(nlsfitSS)$coefficients[3]

    # 2. Derive estimates from the model
    # ---------------------------------
    # current fraction of saturation achieved
    total_loops <- loop_counts[is.na(Replicate), N]
    cur_sat_frac <- total_loops / Asym

    # get model-predicted loop and sample numbers at various levels of saturation
    saturation_ests <- data.table(
        Frac_Saturation = c(0.5, 0.9, 0.95, 0.99)
    )
    saturation_ests[, N_Loops := Asym * Frac_Saturation]
    saturation_ests[, Seq_Depth := -log((Frac_Saturation - 1) * Asym / (R0 - Asym)) / exp(lrc) ]

    # get model-predicted loop and saturation numbers based on the number samples
    depth_ests <- data.table(
        Seq_Depth = seq(50000000, max_depth, 50000000)
    )
    depth_ests[, Est_N_Loops := pmax(
        Asym + (R0 - Asym) * exp(-exp(lrc) * Seq_Depth),
        0
    )]
    depth_ests[, Frac_Saturation := Est_N_Loops / Asym]

    # return all objects in a structured list
    return(list(
        "model" = list(
            "Asym" = Asym,
            "R0" = R0,
            "lrc" = lrc,
            "RSS" = RSS,
            "TSS" = TSS,
            "R_Squared" = r2
        ),
        "loops" = list(
            "identified" = total_loops,
            "saturation" = cur_sat_frac
        ),
        "estimates" = list(
            "saturation" = saturation_ests,
            "depth" = depth_ests
        )
    ))
}


# ==============================================================================
# Data
# ==============================================================================
loginfo("Loading data")
# load metadata
metadata <- fread("config.tsv", sep = "\t", header = TRUE)
metadata <- metadata[Include == "Yes", .SD]
SAMPLES <- list(
    "all" = metadata[, SampleID],
    "benign" = metadata[Type == "Benign", SampleID],
    "tumour" = metadata[Type == "Malignant", SampleID],
    "primary" = metadata[Source == "Primary", SampleID],
    "T2E" = metadata[Type == "Malignant" & T2E == "Yes", SampleID],
    "NonT2E" = metadata[Type == "Malignant" & T2E == "No", SampleID]
)

# load loop calls
sim_loops <- rbindlist(lapply(
    SAMPLES[["primary"]],
    function(s) {
        dt1 <- fread(
            file.path(
                LOOP_DIR,
                paste0(s, ".all-iterations.loops.tsv")
            ),
            sep = "\t",
            header = TRUE
        )
        dt1[, SampleID := s]
        return(dt1)
    }
))

sim_loops[, `:=`(
    SampleID = factor(SampleID, levels = SAMPLES[["primary"]]),
    Seq_Depth = factor(Seq_Depth, levels = DEPTHS, ordered = TRUE),
    Replicate = factor(Replicate, levels = 1:10, ordered = TRUE)
)]

# load full-depth loop calls
full_loops <- rbindlist(lapply(
    SAMPLES[["primary"]],
    function(s) {
        dt1 <- fread(
            file.path(
                "..", "..", "results", "2020-10-06_loops", "Loops", "by-sample",
                paste0(s, ".loops.bedpe")
            ),
            sep = "\t",
            header = FALSE,
            col.names = c(
                "chr_x", "start_x", "end_x",
                "chr_y", "start_y", "end_y",
                "anchor_ID_x", "anchor_ID_y", "loop_ID",
                "fdr", "detection_scale"
            )
        )
        dt1[, SampleID := s]
        return(dt1)
    }
))

# merge in sequencing depth information that's stored in the metadata
full_loops <- merge(
    x = full_loops,
    y = metadata[, .SD, .SDcols = c("SampleID", "Seq_Depth")],
    by = "SampleID"
)


# ==============================================================================
# Analysis
# ==============================================================================
loginfo("Counting interaction calls")
# set the keys for a cross-join operations that includes all factors,
# even the ones without loops
setkeyv(sim_loops, c("SampleID", "Seq_Depth", "Replicate"))
# count the loops in a given (sample, depth, replicate) tuple
sim_loop_counts <- sim_loops[,
    .N,
    by = c("SampleID", "Seq_Depth", "Replicate")
][
    # ensure that tuples with N = 0 are still counted
    CJ(levels(SampleID), levels(Seq_Depth), levels(Replicate))
]
# previous CJ(...) command returns NA when the tuple isn't found
# switch these to 0 for proper counting
sim_loop_counts[is.na(N), N := 0]

# marking the full depth samples as Replicate = NA for easy identification in
# saturation analysis
full_loop_counts <- full_loops[,
    .(
        Replicate = NA,
        N = .N
    ),
    by = c("SampleID", "Seq_Depth")
]

# count all interactions called in each sample
all_loop_counts <- rbindlist(list(
    sim_loop_counts,
    full_loop_counts
))

# set Seq_Depth to a numeric value instead of a factor
# (if this isn't done, there's a problem with the SSasymp function)
all_loop_counts[, Seq_Depth := as.integer(Seq_Depth)]

# perform the saturation analyses
# exclude Benign-Prostate-1664855 from this saturation analysis
# for some reason (that I have yet to figure out), there is a singular gradient
# for this sample when fitting the nls model.
# To avoid this issue for now, just use the other 16 samples.
nls_samples <- setdiff(SAMPLES[["primary"]], "Benign-Prostate-1664855")
saturation_models <- lapply(
    nls_samples,
    function(s) {
        estimate_saturation(all_loop_counts[SampleID == s])
    }
)
names(saturation_models) <- nls_samples

# combine tables to be saved into a major table for each model
depth_ests <- rbindlist(lapply(
    names(saturation_models),
    function(model) {
        dt <- saturation_models[[model]]$estimates$depth
        dt[, Seq_Depth := as.integer(Seq_Depth)]
        dt[, Model := model]
        return(dt)
    }
))
saturation_ests <- rbindlist(lapply(
    names(saturation_models),
    function(model) {
        dt <- saturation_models[[model]]$estimates$saturation
        dt[, Seq_Depth := as.integer(Seq_Depth)]
        dt[, Model := model]
        return(dt)
    }
))


# ==============================================================================
# Save tables
# ==============================================================================
loginfo("Saving tables")
fwrite(
    depth_ests,
    file.path(
        SAT_DIR,
        "loop-saturation.model-estimates.tsv"
    ),
    sep = "\t",
    col.names = TRUE
)
fwrite(
    saturation_ests,
    file.path(
        SAT_DIR,
        "loop-saturation.saturation-estimates.tsv"
    ),
    sep = "\t",
    col.names = TRUE
)
fwrite(
    all_loop_counts,
    file.path(
        SAT_DIR,
        "loop-saturation.iterations.tsv"
    ),
    sep = "\t",
    col.names = TRUE
)

# ==============================================================================
# Plots
# ==============================================================================
loginfo("Generating plots")
for (model in names(saturation_models)) {
    gg <- (
        ggplot()
        + geom_point(
            data = all_loop_counts[SampleID == model],
            mapping = aes(x = Seq_Depth, y = N, colour = Seq_Depth),
            position = position_jitter(width = 0.2, height = 0),
            alpha = 0.1
        )
        + geom_boxplot(
            data = all_loop_counts[SampleID == model],
            mapping = aes(x = Seq_Depth, y = N, group = Seq_Depth, fill = Seq_Depth),
            outlier.shape = NA,
            alpha = 0.5
        )
        # asymptotic curve fit
        + geom_path(
            data = depth_ests[Model == model],
            aes(x = Seq_Depth, y = Est_N_Loops),
            linetype = "dashed"
        )
        # label for model
        + geom_text(
            aes(
                x = saturation_models[[ model ]]$estimates$saturation[Frac_Saturation >= 0.95, mean(Seq_Depth)],
                y = 0.9 * saturation_models[[ model ]]$model$Asym,
                label = "Model Fit"
            ),
            vjust = 0,
            hjust = 0.5
        )
        # add asymptotic value
        + geom_hline(
            aes(yintercept = c(
                saturation_models[[ model ]]$model$Asym,
                saturation_models[[ model ]]$loops$identified
            )),
            linetype = "dashed",
            colour = "#1e90ff"

        )
        # labels for horizontal lines
        + geom_text(
            aes(
                x = c(0, 0),
                y = c(
                    saturation_models[[ model ]]$loops$identified,
                    saturation_models[[ model ]]$model$Asym
                ),
                label = c(
                    paste0(
                        "Interactions Detected: ",
                        saturation_models[[ model ]]$loops$identified,
                        " (", 100 * round(saturation_models[[ model ]]$loops$saturation, 3),
                        "%)"
                    ),
                    paste("Estimated Number of Interactions:", round(saturation_models[[ model ]]$model$Asym))
                )
            ),
            vjust = -1,
            hjust = 0,
            colour = "#1e90ff"
        )
        # add bars for number of samples required to reach saturation
        + geom_vline(
            aes(xintercept = saturation_models[[ model ]]$estimates$saturation[, round(Seq_Depth)]),
            linetype = "dashed",
            colour = "#66cdaa"
        )
        # labels for vertical lines
        + geom_text(
            data = saturation_models[[ model ]]$estimates$saturation,
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
                saturation_models[[ model ]]$estimates$saturation[
                    Frac_Saturation == 0.99,
                    round(Seq_Depth) + 1
                ]
            ),
            breaks = seq(500000000, 2e9, 500000000),
            labels = seq(500, 2e3, 500)
        )
        + scale_y_continuous(
            name = "Number of interactions",
            limits = c(0, 1.05 * saturation_models[[ model ]]$model$Asym)
        )
        + guides(fill = FALSE, colour = FALSE)
        + theme_minimal()
    )
    savefig(gg, file.path(PLOT_DIR, paste0("loop-saturation.", model)))
}

gg <- (
    ggplot()
    # asymptotic curve fit
    + geom_path(
        data = depth_ests,
        aes(x = Seq_Depth, y = Est_N_Loops, colour = Model),
        linetype = "dashed"
    )
    # # add asymptotic value
    # + geom_hline(
    #     aes(yintercept = c(
    #         saturation_models[[ model ]]$model$Asym,
    #         saturation_models[[ model ]]$loops$identified
    #     )),
    #     linetype = "dashed",
    #     colour = "#1e90ff"
    # )
    + scale_x_continuous(
        name = "Sequencing Depth"
    )
    + scale_y_continuous(
        name = "Number of interactions"
    )
    # + guides(fill = FALSE, colour = FALSE)
    + theme_minimal()
)
savefig(gg, file.path(PLOT_DIR, "loop-saturation.all"))
