# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("Matrix"))
suppressMessages(library("ggplot2"))
suppressMessages(library("dryhic"))
suppressMessages(library("dplyr"))
suppressMessages(library("argparse"))

if (!interactive()) {
    PARSER <- argparse::ArgumentParser(
        description = "Plot Hi-C contact maps"
    )
    PARSER$add_argument(
        "raw",
        type = "character",
        help = "Sparse matrix file with raw read counts"
    )
    PARSER$add_argument(
        "norm",
        type = "character",
        help = "Sparse matrix file with normalized read counts"
    )
    PARSER$add_argument(
        "bin_signal",
        type = "character",
        help = "TSV file with normalized bin signal information"
    )
    PARSER$add_argument(
        "domains",
        type = "character",
        help = "TSV file with called domain information"
    )
    PARSER$add_argument(
        "res",
        type = "integer",
        help = "Bin size resolution"
    )
    PARSER$add_argument(
        "-p", "--prefix",
        type = "character",
        help = "Output file prefix",
        default = "contacts"
    )
    PARSER$add_argument(
        "-c", "--chr",
        type = "character",
        help = "Specific chromosome to plot. Otherwise, a plot will be made for each."
    )
    ARGS <- PARSER$parse_args()
}

CHRS = paste0("chr", c(1:22, "X", "Y"))

# update prefix to include resolution
ARGS$prefix = paste0(ARGS$prefix, ".", as.integer(ARGS$res), "bp")


# ==============================================================================
# Data
# ==============================================================================
bins = fread(ARGS$bin_signal, sep = "\t", header = TRUE)
domains = fread(ARGS$domains, sep = "\t", header = TRUE)
mtx_raw = readMM(ARGS$raw)
mtx_norm = readMM(ARGS$norm)

bin_names = bins[, paste0(chr, ':', from.coord)]
dimnames(mtx_raw) = list(bin_names, bin_names)
dimnames(mtx_norm) = list(bin_names, bin_names)


# ==============================================================================
# Plots
# ==============================================================================
if (!is.null(ARGS$chr)) {
    chrs_to_plot = ARGs$chr
} else {
    chrs_to_plot = CHRS
}

cat("Plotting\n")
for (plot_chr in chrs_to_plot) {
    cat("  ", plot_chr, "\n")
    bin_idx = bins[, which(chr == plot_chr)]
    coord_lims = c(
        bins[bin_idx, ][1, from.coord],
        bins[bin_idx, ][length(bin_idx), from.coord] + ARGS$res
    )
    domains_tuples = domains[chr == plot_chr & tag == "domain", .(from.id, to.id)]
    pdf(
        paste(ARGS$prefix, plot_chr, "raw", "pdf", sep = "."),
        height = 12,
        width = 12
    )
    plot_matrix(
        mat = as.matrix(mtx_raw[bin_idx, bin_idx]),
        coord = coord_lims,
        resolution = ARGS$res
    )
    dev.off()

    pdf(
        paste(ARGS$prefix, plot_chr, "ice-corrected", "pdf", sep = "."),
        height = 12,
        width = 12
    )
    plot_matrix(
        mat = as.matrix(mtx_norm[bin_idx, bin_idx]),
        coord = coord_lims,
        resolution = ARGS$res
    )
    dev.off()

    pdf(
        paste(ARGS$prefix, plot_chr, "ice-corrected-with-TADs", "pdf", sep = "."),
        height = 12,
        width = 12
    )
    plot_matrix(
        mat = as.matrix(mtx_norm[bin_idx, bin_idx]),
        coord = coord_lims,
        resolution = ARGS$res,
        tads = domains_tuples
    )
    dev.off()
}


warnings()