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
PLOT_DIR <- file.path(RES_DIR, "Plots")

# ==============================================================================
# Functions
# ==============================================================================
#' Use bootstraps and an exponential fit to estimate the total number of loops
#' present in a batch of samples
#'
#' @param sample_loops Loops (or IDs for the loops) present in each sample
#' @param n_iters Number of bootstraps to generate
#' @param max_samples maximum number of samples to estimate the model out to
#' @return structured list
estimate_saturation <- function(sample_loops, n_iters = 100, max_samples = 200) {
    # 1. Generate bootstraps
    # ---------------------------------
    sample_names <- names(sample_loops)
    n_loops <- lengths(sample_loops)
    # number of samples
    n_samples = length(sample_loops)
    # create the results list
    data_iter = matrix(ncol = n_iters, nrow = n_samples)

    for (j in 1:n_iters) {
        # randomly order the samples
        reordered_loops <- sample_loops[sample(1:n_samples)]

        # create cumulative set of loop IDs
        reordered_loop_IDs <- vector("list", length = n_samples)
        reordered_loop_IDs[[1]] <- reordered_loops[[1]]
        # find union of loop IDs from the combined set of i samples
        # doing this in the for loop has a "cumulative union" effect
        for (i in 2:n_samples) {
            reordered_loop_IDs[[i]] <- union(
                reordered_loop_IDs[[i - 1]],
                reordered_loops[[i]]
            )
        }
        # count how many there are in this iteration
        reordered_loop_ID_counts <- lengths(reordered_loop_IDs)

        # store data in saturation curve matrix as the j-th column
        data_iter[, j] <- reordered_loop_ID_counts
    }

    # 2. Calculate summary statistics
    # ---------------------------------
    bootstrap_summary <- data.table(
        N_Samples = 1:n_samples,
        N_Loops_Mean = rowMeans(data_iter),
        N_Loops_SEM = rowSds(data_iter) / sqrt(n_iters)
    )

    # coerce into a data.table for plotting the bootstrap results themselves
    bootstraps <- as.data.table(cbind(1:n_samples, data_iter))
    colnames(bootstraps) <- c("N_Samples", 1:n_iters)
    bootstraps_long <- melt(
        bootstraps,
        id.vars = "N_Samples",
        variable.name = "Iteration",
        value.name = "N_Loops"
    )
    bootstraps_long[, Iteration := as.numeric(Iteration)]

    # 3. Fit exponential model to the bootstraps
    # ---------------------------------
    # use self starting model for asymptotic data SSasymp
    # it follows the formula: y = Asym + (R0 - Asym) * exp(-exp(lrc) * x) where
    # x     a numeric vector of values at which to evaluate the model
    # Asym  a numeric parameter for the upper asymptotic value of the model.
    # R0    a numeric parameter for the response when input is zero
    # lrc	a numeric parameter for the natural log of the rate constant
    nlsfitSS <- nls(
        N_Loops_Mean ~ SSasymp(N_Samples, Asym, R0, lrc),
        data = bootstrap_summary,
        control = list(maxiter = 10000)
    )

    # predict the future values with increasing numbers of samples
    pred <- predict(nlsfitSS, list(N_Samples = seq(1, max_samples)))

    # check that the model is a good fit for the data
    #   Residual sum of squares
    RSS <- sum(residuals(nlsfitSS)^2)
    #   Total sum of squares
    TSS <- bootstrap_summary[, sum((N_Loops_Mean - mean(N_Loops_Mean))^2)]
    #   R-squared measure (this should be close to 1)
    r2 <- 1 - (RSS/TSS)

    # extract key parameters from the model
    Asym = summary(nlsfitSS)$coefficients[1]
    R0 = summary(nlsfitSS)$coefficients[2]
    lrc = summary(nlsfitSS)$coefficients[3]

    # 4. Derive estimates from the model
    # ---------------------------------
    # current fraction of saturation achieved
    total_loops <- bootstrap_summary[N_Samples == n_samples, N_Loops_Mean]
    cur_sat_frac <- total_loops / Asym

    # get model-predicted loop and sample numbers at various levels of saturation
    saturation_ests <- data.table(
        Frac_Saturation = c(0.5, 0.9, 0.95, 0.99)
    )
    saturation_ests[, N_Loops := Asym * Frac_Saturation]
    saturation_ests[, N_Samples := -log((Frac_Saturation - 1) * Asym / (R0 - Asym)) / exp(lrc) ]

    # get model-predicted loop and saturation numbers based on the number samples
    sample_ests <- data.table(
        N_Samples = 1:max_samples
    )
    sample_ests[, Est_N_Loops := Asym + (R0 - Asym) * exp(-exp(lrc) * N_Samples)]
    sample_ests[, Frac_Saturation := Est_N_Loops / Asym]

    # return all objects in a structured list
    return(list(
        "bootstraps" = bootstraps_long,
        "bootstrap_summary" = bootstrap_summary,
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
            "samples" = sample_ests
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
    "T2E" = metadata[Type == "Malignant" & T2E == "Yes", SampleID],
    "NonT2E" = metadata[Type == "Malignant" & T2E == "No", SampleID]
)

# load loop calls
loops <- fread(
    file.path(LOOP_DIR, "merged-loops.tsv"),
    sep = "\t",
    header = TRUE
)


# ==============================================================================
# Analysis
# ==============================================================================
loginfo("Generating bootstraps and fitting models")
# create peaks list
loops_per_sample <- lapply(
    SAMPLES$all,
    function(s) {
        loops[SampleID == s, loop_ID]
    }
)
names(loops_per_sample) <- SAMPLES$all

loops_tumour <- loops_per_sample[SAMPLES[["tumour"]]]
loops_benign <- loops_per_sample[SAMPLES[["benign"]]]
loops_t2e <- loops_per_sample[SAMPLES[["T2E"]]]
loops_nont2e <- loops_per_sample[SAMPLES[["NonT2E"]]]

# perform the saturation analyses
saturation_models <- list(
    "All" = estimate_saturation(loops_per_sample),
    "Tumour_Only" = estimate_saturation(loops_tumour),
    "Benign_Only" = estimate_saturation(loops_benign),
    "T2E_Tumour_Only" = estimate_saturation(loops_t2e),
    "NonT2E_Tumour_Only" = estimate_saturation(loops_nont2e)
)

# combine tables to be saved into a major table for each model
sample_ests <- rbindlist(lapply(
    names(saturation_models),
    function(model) {
        dt <- saturation_models[[model]]$estimates$samples
        dt[, Model := model]
        return(dt)
    }
))
saturation_ests <- rbindlist(lapply(
    names(saturation_models),
    function(model) {
        dt <- saturation_models[[model]]$estimates$saturation
        dt[, Model := model]
        return(dt)
    }
))
bootstrap_summaries <- rbindlist(lapply(
    names(saturation_models),
    function(model) {
        dt <- saturation_models[[model]]$bootstrap_summary
        dt[, Model := model]
        return(dt)
    }
))
bootstraps_long <- rbindlist(lapply(
    names(saturation_models),
    function(model) {
        dt <- saturation_models[[model]]$bootstraps
        dt[, Model := model]
        return(dt)
    }
))

# ==============================================================================
# Save tables
# ==============================================================================
loginfo("Saving tables")
fwrite(
    sample_ests,
    file.path(
        LOOP_DIR,
        "loop-saturation.model-estimates.tsv"
    ),
    sep = "\t",
    col.names = TRUE
)
fwrite(
    saturation_ests,
    file.path(
        LOOP_DIR,
        "loop-saturation.saturation-estimates.tsv"
    ),
    sep = "\t",
    col.names = TRUE
)
fwrite(
    bootstrap_summaries,
    file.path(
        LOOP_DIR,
        "loop-saturation.bootstrap-summary.tsv"
    ),
    sep = "\t",
    col.names = TRUE
)
fwrite(
    bootstraps_long,
    file.path(
        LOOP_DIR,
        "loop-saturation.bootstraps.tsv"
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
            data = bootstraps_long[Model == model],
            mapping = aes(x = N_Samples, y = N_Loops, colour = N_Samples),
            position = position_jitter(width = 0.2, height = 0),
            alpha = 0.1
        )
        + geom_boxplot(
            data = bootstraps_long[Model == model],
            mapping = aes(x = N_Samples, y = N_Loops, group = N_Samples, fill = N_Samples),
            outlier.shape = NA,
            alpha = 0.5
        )
        # asymptotic curve fit
        + geom_path(
            data = sample_ests[Model == model],
            aes(x = N_Samples, y = Est_N_Loops),
            linetype = "dashed"
        )
        # label for model
        + geom_text(
            aes(
                x = saturation_models[[ model ]]$estimates$saturation[Frac_Saturation >= 0.95, mean(N_Samples)],
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
                        "Loops Detected: ",
                        saturation_models[[ model ]]$loops$identified,
                        " (", 100 * round(saturation_models[[ model ]]$loops$saturation, 3),
                        "%)"
                    ),
                    paste("Estimated Number of Loops:", round(saturation_models[[ model ]]$model$Asym))
                )
            ),
            vjust = -1,
            hjust = 0,
            colour = "#1e90ff"
        )
        # add bars for number of samples required to reach saturation
        + geom_vline(
            aes(xintercept = saturation_models[[ model ]]$estimates$saturation[, round(N_Samples)]),
            linetype = "dashed",
            colour = "#66cdaa"
        )
        # labels for vertical lines
        + geom_text(
            data = saturation_models[[ model ]]$estimates$saturation,
            aes(
                x = round(N_Samples),
                y = 0,
                label = paste0(
                    100 * Frac_Saturation,
                    "%: ",
                    round(N_Samples),
                    " samples"
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
            name = "Number of samples",
            limits = c(
                0,
                # estimated number of samples required to reach 99% saturation of loop calls
                saturation_models[[ model ]]$estimates$saturation[
                    Frac_Saturation == 0.99,
                    round(N_Samples) + 1
                ]
            )
        )
        + scale_y_continuous(
            name = "Number of loops",
            limits = c(0, 1.05 * saturation_models[[ model ]]$model$Asym)
        )
        + guides(fill = FALSE, colour = FALSE)
        + theme_minimal()
    )
    savefig(gg, file.path(PLOT_DIR, paste0("loop-saturation.", model)))
}
