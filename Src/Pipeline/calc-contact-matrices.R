# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("dryhic"))
suppressMessages(library("Matrix"))
suppressMessages(library("argparse"))

if (!interactive()) {
    PARSER <- argparse::ArgumentParser(
        description = "Generate raw and normalized contact matrices from a BAM at a desired resolution"
    )
    PARSER$add_argument(
        "bam",
        type = "character",
        help = "Input filtered BAM"
    )
    PARSER$add_argument(
        "res",
        type = "integer",
        help = "Resolution (size) of bins (bp)"
    )
    PARSER$add_argument(
        "-p", "--prefix",
        type = "character",
        help = "Prefix for output files",
        default = "contacts"
    )
    ARGS <- PARSER$parse_args()
} else {
    ARGS = list(
        bam = "Aligned/PCa13848_R1_2.hicup.name-sorted.bam",
        res = 50000,
        prefix = "Contacts/PCa13848"
    )
}

CHRS = paste0("chr", c(1:22, "X", "Y"))
# update prefix to include resolution
ARGS$prefix = paste0(ARGS$prefix, ".", ARGS$res, "bp")


# ==============================================================================
# Data
# ==============================================================================
# bins generated have 0-indexed start, 1-indexed end
bins = make_bins(inbam = ARGS$bam, resolution = ARGS$res)
# BAMs have 1-indexed start
mtx_raw = get_contacts_matrix(
    inbam = ARGS$bam,
    resolution = ARGS$res,
    pos = bins$bin,
    filtin = 3, # paired and each mate is aligned
    filtex = 12 # unmapped or mate is unmapped
)

mtx_ice = ICE(mtx_raw, 30, verbose = FALSE)

bins_bed = data.table(
    chr = bins$chr,
    from.coord = as.integer(bins$pos),
    to.coord = as.integer(bins$pos + ARGS$res),
    name = bins$bin
)

# ==============================================================================
# Analysis
# ==============================================================================
# call TADs using TopDom
tadcalls = list(
    binSignal = NULL,
    domain = NULL,
    bed = NULL
)
for (chr in CHRS) {
    cat(chr, "\n")
    # get all bins on current chromosome
    bin_idx = which(bins$chr == chr)
    # convert sparse matrix to dense for processing
    mtx_ice_dense = as.matrix(mtx_ice[bin_idx, bin_idx])
    # use TopDom
    calls = TopDom(
        matrix.data = mtx_ice_dense,
        bins = bins_bed[bin_idx],
        window.size = 5
    )
    # concatenate results
    tadcalls$binSignal = rbind(tadcalls$binSignal, calls$binSignal)
    tadcalls$domain = rbind(tadcalls$domain, calls$domain)
    tadcalls$bed = rbind(tadcalls$bed, calls$bed)
}

# convert to data.tables for saving
domains = as.data.table(tadcalls$domain)
domains_bed = as.data.table(tadcalls$bed)

# ==============================================================================
# Save Data
# ==============================================================================
writeMM(
    mtx_raw,
    paste(ARGS$prefix, "raw", "sparse", "mtx", sep = ".")
)
writeMM(
    mtx_ice,
    paste(ARGS$prefix, "ice-corrected", "sparse", "mtx", sep = ".")
)

fwrite(
    bins_bed,
    paste(ARGS$prefix, "bins", "bed", sep = "."),
    sep = "\t",
    col.names = FALSE
)

fwrite(
    tadcalls$binSignal,
    paste(ARGS$prefix, "bins-signal", "tsv", sep = "."),
    sep = "\t",
    col.names = TRUE
)

fwrite(
    domains,
    paste(ARGS$prefix, "domains", "tsv", sep = "."),
    sep = "\t",
    col.names = TRUE
)

fwrite(
    domains_bed,
    paste(ARGS$prefix, "domains", "bed", sep = "."),
    sep = "\t",
    col.names = FALSE
)

# ==============================================================================
# Plots
# ==============================================================================

bin_idx = which(bins$chr == "chr21")
coord_lims = c(
    bins[bin_idx, ][1, "pos"],
    bins[bin_idx, ][length(bin_idx), "pos"] + ARGS$res
)

domains_tuples = domains[chr == "chr21" & tag == "domain", .(from.id, to.id)]

pdf(
    paste(ARGS$prefix, "raw", "pdf", sep = "."),
    height = 12,
    width = 12
)
plot_matrix(
    mat = mtx_raw[which(bins$chr == "chr21"), which(bins$chr == "chr21")],
    coord = coord_lims,
    resolution = ARGS$res
)
dev.off()

pdf(
    paste(ARGS$prefix, "ice-corrected", "pdf", sep = "."),
    height = 12,
    width = 12
)
plot_matrix(
    mat = mtx_ice[which(bins$chr == "chr21"), which(bins$chr == "chr21")],
    coord = coord_lims,
    resolution = ARGS$res
)
dev.off()

pdf(
    paste(ARGS$prefix, "ice-corrected-with-TADs", "pdf", sep = "."),
    height = 12,
    width = 12
)
plot_matrix(
    mat = mtx_ice[which(bins$chr == "chr21"), which(bins$chr == "chr21")],
    coord = coord_lims,
    resolution = ARGS$res,
    tads = domains_tuples
)
dev.off()
