# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("DiffBind"))
suppressMessages(library("ggplot2"))
source("../2020-02-19_chromoplexy/plotting-helper.R")

BAM_DIR <- file.path("..", "..", "Data", "Processed", "2019-05-03_PCa-H3K27ac-peaks", "BAMs")
PEAK_DIR <- file.path("..", "..", "Data", "Processed", "2019-05-03_PCa-H3K27ac-peaks", "Peaks")
TEST_DIR <- file.path("Acetylation", "T2E")
PLOT_DIR <- file.path("Plots", "T2E")


# ==============================================================================
# Data
# ==============================================================================
# load sample metadata
metadata <- fread(
    "config.tsv",
    sep = "\t",
    header = TRUE
)
low_qual_samples <- metadata[Include == "No"]
metadata <- metadata[Include == "Yes"]
SAMPLES <- metadata$SampleID

# load DiffBind data
dba_comp <- readRDS(file.path(TEST_DIR, "dba_comp.rds"))
local_comp <- fread(file.path(TEST_DIR, "t2e.local.tsv"))

# ==============================================================================
# Plots
# ==============================================================================
gg_conc <- (
    ggplot(data = local_comp[order(start, end)])
    + geom_smooth(aes(x = start, y = Conc_T2E), method = "loess", colour = "#0000cd")
    + geom_smooth(aes(x = start, y = Conc_NonT2E), method = "loess", colour = "#ff8c00")
    + labs(x = "chr21 position", y = "H3K27ac Concentration")
    + theme_minimal()
)
savefig(gg_conc, file.path(PLOT_DIR, "position"))
