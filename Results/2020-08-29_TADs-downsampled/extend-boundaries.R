# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("logging"))

loginfo("Loading packages")
suppressMessages(library("data.table"))
suppressMessages(library("argparse"))

if (!interactive()) {
    PARSER <- argparse::ArgumentParser(
        description = "Call TADs on a contact matrix using TopDom"
    )
    PARSER$add_argument(
        "boundaries",
        type = "character",
        help = "Path to TAD boundary calls"
    )
    PARSER$add_argument(
        "extension",
        type = "integer",
        help = "Number of basepairs to extend the boundary on either side"
    )
    PARSER$add_argument(
        "output",
        type = "character",
        help = "Output file to save to"
    )
    ARGS <- PARSER$parse_args()
} else {
    ARGS <- list(
        boundaries = "Aggregated-TADs/PCa13266.300000000.res_40000bp.agg-boundaries.tsv",
        extension = 200000,
        output = "CTCF/Extended/PCa13266.300000000.res_40000bp.agg-boundaries.extended-200000bp.bed"
    )
}


# ==============================================================================
# Data
# ==============================================================================
loginfo("Loading Data\n")

# load aggregated boundary calls from each sample
boundaries <- fread(ARGS$boundaries, header = TRUE, sep = "\t")

# ==============================================================================
# Analysis
# ==============================================================================
loginfo("Extending boundaries\n")
ext_boundaries <- boundaries[
    ,
    .(
        "chr_ext" = chr,
        # force integer to ensure no scientific notation used in text file
        "start_ext" = as.integer(pmax(0, pos - ARGS$extension)),
        "end_ext" = as.integer(pos + ARGS$extension),
        "chr_bound" = chr,
        "pos_bound" = as.integer(pos)
    )
]

# ==============================================================================
# Save data
# ==============================================================================
fwrite(
    ext_boundaries,
    ARGS$output,
    sep = "\t",
    col.names = FALSE
)