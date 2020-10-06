# ==============================================================================
# Meta
# ==============================================================================
# plot-insulation
# --------------------------------------
# Description: Make plots related to chromatin insulation
# Author: James Hawley


# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
suppressMessages(library("logging"))

# ==============================================================================
# Functions
# ==============================================================================
#' Save figures in multiple formats
#'
#' @param gg ggplot object
#' @param prefix Prefix for output file
#' @param ext Output extensions
#' @param dpi DPI resolution
savefig <- function(gg, prefix, ext = c("png", "pdf"), width = 20, height = 12, dpi = 400) {
    for (e in ext) {
        ggsave(
            paste(prefix, e, sep = "."),
            gg,
            height = height,
            width = width,
            units = "cm",
            dpi = dpi
        )
    }
}

# ==============================================================================
# Data
# ==============================================================================
loginfo("Loading Data")
# load metadata
metadata <- fread("config.tsv", sep = "\t", header = TRUE)
metadata <- metadata[Include == "Yes"]
SAMPLES <- metadata[, SampleID]
TUMOUR_SAMPLES <- metadata[Source == "Primary" & Type == "Malignant", SampleID]
BENIGN_SAMPLES <- metadata[Source == "Primary" & Type == "Benign", SampleID]
PRIMARY_SAMPLES <- metadata[Source == "Primary", SampleID]
n_tumour <- length(TUMOUR_SAMPLES)
n_benign <- length(BENIGN_SAMPLES)

# load TADs
insulation <- rbindlist(lapply(
    PRIMARY_SAMPLES,
    function(s) {
        dt <- fread(
            file.path("Insulation", paste0(s, ".300000000.insulation.tsv")),
            sep = "\t",
            header = TRUE
        )
        # remove "bad bins"
        dt <- dt[is_bad_bin == FALSE, .SD]
        dt[, SampleID := s]
        return(dt)
    }
))

# merge metadata with insulation scores
insulation <- merge(
    x = insulation,
    y = metadata[, .SD, .SDcols = c("SampleID", "Type")],
    by = "SampleID"
)



# ==============================================================================
# Analysis
# ==============================================================================
# calculate median insulation score
med_insulation <- insulation[,
    .(log2_insulation = median(log2_insulation_score_120000, na.rm = TRUE)),
    keyby = c("chrom", "start", "end", "Type")
]

# convert to wide format for plotting
wide_ins <- dcast(
    med_insulation,
    chrom + start + end ~ Type,
    value.var = "log2_insulation"
)

rho <- wide_ins[, cor(Benign, Malignant, method = "spearman", use = "pairwise.complete.obs")]

# ==============================================================================
# Plots
# ==============================================================================
loginfo("Plotting data")
gg <- (
    ggplot(data = wide_ins)
    + geom_point(aes(x = Benign, y = Malignant), alpha = 0.1)
    + labs(
        x = expression("Median Benign Insulation (" * log[2] * ")"),
        y = expression("Median Tumour Insulation (" * log[2] * ")")
    )
    + annotate(
        geom = "text",
        x = 0.5,
        y = -2,
        label = bquote("Spearman's " * rho * " = " * .(rho))
    )
    + scale_fill_viridis_c()
    + theme_minimal()
)
savefig(gg, "Plots/primary/insulation-correlation", height = 10, width = 10)
