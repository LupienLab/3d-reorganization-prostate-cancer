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
        "-l", "--locus",
        type = "character",
        help = "Specific genomic locus to plot. Format {chr}:{start}-{end}. If not specificed, a plot will be made for each chromosome, separately."
    )
    ARGS <- PARSER$parse_args()
} else {
    ARGS = list(
        raw = "PCa13848.50000bp.raw.sparse.mtx",
        norm = "PCa13848.50000bp.ice-corrected.sparse.mtx",
        bin_signal = "PCa13848.50000bp.bins-signal.tsv",
        domains = "PCa13848.50000bp.domains.tsv",
        res = 50000,
        locus = "chr21:38350000-41550000"
    )
}

# check that locus format meets specification
if (!is.null(ARGS$locus)) {
    locus_split = unlist(strsplit(ARGS$locus, "[:-]"))
    plot_locus = data.table(
        chr = locus_split[1],
        start = as.integer(locus_split[2]),
        end = as.integer(locus_split[3])
    )
} else {
    plot_locus == NULL
}

CHRS = data.table(
    chr = paste0("chr", c(1:22, "X", "Y")),
    start = 0,
    end = 1e9
)

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
if (!is.null(plot_locus)) {
    locus_to_plot = plot_locus
} else {
    locus_to_plot = CHRS
}

cat("Plotting\n")
for (i in 1:dim(locus_to_plot)[1]) {
    plot_chr = locus_to_plot[i, chr]
    plot_start = locus_to_plot[i, start]
    plot_end = locus_to_plot[i, end]
    cat("  ", paste0(plot_chr, ":", plot_start, "-", plot_end), "\n")
    bin_idx = bins[, which(chr == plot_chr & from.coord >= plot_start & to.coord <= plot_end)]
    coord_lims = c(
        bins[bin_idx, ][1, from.coord],
        bins[bin_idx, ][length(bin_idx), from.coord] + ARGS$res
    )
    domains_tuples = domains[chr == plot_chr & from.coord >= plot_start & from.coord <= plot_end & tag == "domain", .(from.id, to.id)]
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