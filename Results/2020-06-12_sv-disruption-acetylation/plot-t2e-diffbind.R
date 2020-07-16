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

# load GENCODE data
gencode <- fread(
    file.path("..", "..", "Data", "External", "GENCODE", "gencode.v33.genes.sorted.bed"),
    sep = "\t",
    header = FALSE,
    col.names = c("chr", "start", "end", "name", "score", "strand", "gene_id")
)

# ==============================================================================
# Plots
# ==============================================================================
local_comp_long <- melt(
    local_comp,
    id.vars = c("chr", "start", "end"),
    measure.vars = c("Conc_T2E", "Conc_NonT2E"),
    variable.name = "T2E",
    value.name = "H3K27ac"
)

gg_conc <- (
    ggplot()
    + geom_path(
        data = local_comp_long[order(start, end)],
        mapping = aes(x = start, y = H3K27ac, colour = T2E)
    )
    + geom_segment(
        data = gencode[name %in% c("ERG", "TMPRSS2")],
        mapping = aes(x = start, xend = end, y = 0, yend = 0),
        colour = "#000000"
    )
    + geom_text(
        data = gencode[name %in% c("ERG", "TMPRSS2")],
        mapping = aes(x = (start + end) / 2, y = 0, label = name),
        colour = "#000000",
        hjust = 0.5,
        vjust = -1
    )
    + labs(x = "chr21 position", y = expression(log[2] * "(H3K27ac fold enrichment)"))
    + scale_colour_manual(
        breaks = c("Conc_T2E", "Conc_NonT2E"),
        labels = c("T2E", "Non-T2E"),
        values = c("#0000cd", "#ff8c00"),
        name = "T2E Status"
    )
    + coord_cartesian(
        xlim = c(38e6, 42e6)
    )
    + theme_minimal()
    + theme(
        legend.position = "right"
    )
)
savefig(gg_conc, file.path(PLOT_DIR, "position"))

gg_conc_erg <- gg_conc + coord_cartesian(xlim = c(
    gencode[name == "ERG", start - 1e5],
    gencode[name == "ERG", end + 1e5])
)
savefig(gg_conc_erg, file.path(PLOT_DIR, "position.ERG"))

gg_conc_tmprss2 <- gg_conc + coord_cartesian(xlim = c(
    gencode[name == "TMPRSS2", start - 1e5],
    gencode[name == "TMPRSS2", end + 1e5])
)
savefig(gg_conc_tmprss2, file.path(PLOT_DIR, "position.TMPRSS2"))
