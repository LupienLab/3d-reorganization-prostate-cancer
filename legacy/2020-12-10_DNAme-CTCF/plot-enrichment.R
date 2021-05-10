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
perm_test_all <- readRDS("Overlaps/all-dmrs.permutation-tests.rds")
perm_summ_all <- fread("Overlaps/all-dmrs.permutations.tsv")
lz_all <- readRDS("Overlaps/all-dmrs.local-dependency.rds")
lz_all_dt <- fread("Overlaps/all-dmrs.local-dependency.tsv")

perm_test_up <- readRDS("Overlaps/pca-hyper-dmrs.permutation-tests.rds")
perm_summ_up <- fread("Overlaps/pca-hyper-dmrs.permutations.tsv")
lz_up <- readRDS("Overlaps/pca-hyper-dmrs.local-dependency.rds")
lz_up_dt <- fread("Overlaps/pca-hyper-dmrs.local-dependency.tsv")

perm_test_dn <- readRDS("Overlaps/pca-hypo-dmrs.permutation-tests.rds")
perm_summ_dn <- fread("Overlaps/pca-hypo-dmrs.permutations.tsv")
lz_dn <- readRDS("Overlaps/pca-hypo-dmrs.local-dependency.rds")
lz_dn_dt <- fread("Overlaps/pca-hypo-dmrs.local-dependency.tsv")

CELL_LINES <- names(perm_test_all)

perm_test_data_all <- rbindlist(lapply(
    CELL_LINES,
    function(cl) {
        data.table(
            Cell_Line = cl,
            Iteration = 1:perm_test_all[[cl]]$numOverlaps$ntimes,
            Permuted = perm_test_all[[cl]]$numOverlaps$permuted,
            log2Perm_Obs = log2(perm_test_all[[cl]]$numOverlaps$permuted / perm_test_all[[cl]]$numOverlaps$observed)
        )
    }
))

obs_all <- rbindlist(lapply(
    CELL_LINES,
    function(cl) {
        data.table(
            Cell_Line = cl,
            Observed = perm_test_all[[cl]]$numOverlaps$observed
        )
    }
))

perm_test_data_up <- rbindlist(lapply(
    CELL_LINES,
    function(cl) {
        data.table(
            Cell_Line = cl,
            Iteration = 1:perm_test_up[[cl]]$numOverlaps$ntimes,
            Permuted = perm_test_up[[cl]]$numOverlaps$permuted,
            log2Perm_Obs = log2(perm_test_up[[cl]]$numOverlaps$permuted / perm_test_up[[cl]]$numOverlaps$observed)
        )
    }
))

obs_up <- rbindlist(lapply(
    CELL_LINES,
    function(cl) {
        data.table(
            Cell_Line = cl,
            Observed = perm_test_up[[cl]]$numOverlaps$observed
        )
    }
))


perm_test_data_dn <- rbindlist(lapply(
    CELL_LINES,
    function(cl) {
        data.table(
            Cell_Line = cl,
            Iteration = 1:perm_test_dn[[cl]]$numOverlaps$ntimes,
            Permuted = perm_test_dn[[cl]]$numOverlaps$permuted,
            log2Perm_Obs = log2(perm_test_dn[[cl]]$numOverlaps$permuted / perm_test_dn[[cl]]$numOverlaps$observed)
        )
    }
))

obs_dn <- rbindlist(lapply(
    CELL_LINES,
    function(cl) {
        data.table(
            Cell_Line = cl,
            Observed = perm_test_dn[[cl]]$numOverlaps$observed
        )
    }
))


# ==============================================================================
# Plots
# ==============================================================================
loginfo("Plotting figures")
gg <- (
    ggplot(data = perm_test_data_all)
    + geom_histogram(aes(x = Permuted), bins = 100)
    + geom_vline(
        data = obs_all,
        mapping = aes(xintercept = Observed),
        linetype = "dashed"
    )
    + geom_text(
        data = obs_all,
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
    "Plots/DMR-CTCF.all-dmrs.permutations.hist.png",
    width = 16,
    height = 8,
    units = "cm"
)

gg <- (
    ggplot(data = lz_all_dt)
    + geom_line(aes(x = Shift, y = z))
    + labs(x = "Local shift (bp)", y = "Shifted Z-score")
    + theme_minimal()
    + facet_wrap(~ Cell_Line)
)
ggsave(
    "Plots/DMR-CTCF.all-dmrs.permutations.shifted-z.png",
    height = 12,
    width = 20,
    units = "cm"
)

gg <- (
    ggplot(data = perm_test_data_up)
    + geom_histogram(aes(x = Permuted), bins = 100)
    + geom_vline(
        data = obs_up,
        mapping = aes(xintercept = Observed),
        linetype = "dashed"
    )
    + geom_text(
        data = obs_up,
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
    "Plots/DMR-CTCF.pca-hyper-dmrs.permutations.hist.png",
    width = 16,
    height = 8,
    units = "cm"
)

gg <- (
    ggplot(data = lz_up_dt)
    + geom_line(aes(x = Shift, y = z))
    + labs(x = "Local shift (bp)", y = "Shifted Z-score")
    + theme_minimal()
    + facet_wrap(~ Cell_Line)
)
ggsave(
    "Plots/DMR-CTCF.pca-hyper-dmrs.permutations.shifted-z.png",
    height = 12,
    width = 20,
    units = "cm"
)

gg <- (
    ggplot(data = perm_test_data_dn)
    + geom_histogram(aes(x = Permuted), bins = 100)
    + geom_vline(
        data = obs_dn,
        mapping = aes(xintercept = Observed),
        linetype = "dashed"
    )
    + geom_text(
        data = obs_dn,
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
    "Plots/DMR-CTCF.pca-hypo-dmrs.permutations.hist.png",
    width = 16,
    height = 8,
    units = "cm"
)

gg <- (
    ggplot(data = lz_dn_dt)
    + geom_line(aes(x = Shift, y = z))
    + labs(x = "Local shift (bp)", y = "Shifted Z-score")
    + theme_minimal()
    + facet_wrap(~ Cell_Line)
)
ggsave(
    "Plots/DMR-CTCF.pca-hypo-dmrs.permutations.shifted-z.png",
    height = 12,
    width = 20,
    units = "cm"
)
