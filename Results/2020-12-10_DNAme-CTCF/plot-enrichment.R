# ==============================================================================
# Meta
# ==============================================================================
# plot-enrichment
# --------------------------------------
# Description: Plot enrichment permutations for CTCF binding sites in prostate cancer cell lines
# Author: James Hawley


# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("logging"))

loginfo("Loading packages")
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))

# ==============================================================================
# Functions
# ==============================================================================

# ==============================================================================
# Data
# ==============================================================================
loginfo("Loading data")
perm_test <- readRDS("Overlaps/permutation-tests.rds")
lz <- readRDS("Overlaps/local-dependency.rds")
perm_test_data <- fread("Overlaps/permutations.tsv")

CELL_LINES <- names(perm_test)

perm_tests <- rbindlist(lapply(
    CELL_LINES,
    function(cl) {
        data.table(
            Cell_Line = cl,
            Iteration = 1:perm_test[[cl]]$numOverlaps$ntimes,
            Permuted = perm_test[[cl]]$numOverlaps$permuted,
            log2Perm_Obs = log2(perm_test[[cl]]$numOverlaps$permuted / perm_test[[cl]]$numOverlaps$observed)
        )
    }
))

obs_tests <- rbindlist(lapply(
    CELL_LINES,
    function(cl) {
        data.table(
            Cell_Line = cl,
            Observed = perm_test[[cl]]$numOverlaps$observed
        )
    }
))

lz_dt <- fread("Overlaps/local-dependency.tsv")

# ==============================================================================
# Plots
# ==============================================================================
loginfo("Plotting figures")
gg <- (
    ggplot(data = perm_tests)
    + geom_histogram(aes(x = Permuted), bins = 100)
    + geom_vline(
        data = obs_tests,
        mapping = aes(xintercept = Observed),
        linetype = "dashed"
    )
    + geom_text(
        data = obs_tests,
        mapping = aes(x = Observed, y = 100, label = "Observed"),
        angle = 90,
        hjust = 0.5,
        vjust = -0.5
    )
    + labs(x = "DMR-CTCF Intersections", y = "Frequency")
    + theme_minimal()
    + facet_wrap(~ Cell_Line)
)
ggsave(
    "Plots/DMR-CTCF.permutations.hist.png",
    width = 16,
    height = 8,
    units = "cm"
)

gg <- (
    ggplot(data = lz_dt)
    + geom_line(aes(x = Shift, y = z))
    + labs(x = "Local shift (bp)", y = "Shifted Z-score")
    + theme_minimal()
    + facet_wrap(~ Cell_Line)
)
ggsave(
    "Plots/DMR-CTCF.permutations.shifted-z.png",
    height = 12,
    width = 20,
    units = "cm"
)
