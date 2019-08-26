# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("argparse"))

if (!interactive()) {
    PARSER <- argparse::ArgumentParser(
        description = "Aggregate breakpoint calls from `hic_breakfinder`"
    )
    PARSER$add_argument(
        "-p", "--prefix",
        type = "character",
        help = "Output file prefix. `{prefix}.bedpe` and `{prefix}.tsv` are created",
        default="breakpoints"
    )
    PARSER$add_argument(
        "breakpoints",
        type = "character",
        help = "Breakpoint call files to aggregate",
        nargs="+"
    )
    ARGS <- PARSER$parse_args()
}

# ==============================================================================
# Data
# ==============================================================================
agg_data = rbindlist(lapply(
    ARGS$breakpoints,
    function(bp) {
        # read breakpoint calls
        dt = fread(bp, sep = "\t", header = FALSE, col.names = c("LogOdds", "Chrom1", "Start1", "End1", "Strand1", "Chrom2", "Start2", "End2", "Strand2", "Resolution"))
        # extract sample name from file name
        #   find '.' in file name
        idx = regexpr("\\.", basename(bp))
        #   take everything before first '.'
        samp = substr(basename(bp), 1, idx[1] - 1)
        dt[, Sample := samp]
        return(dt)
    }
))

# ==============================================================================
# Save data
# ==============================================================================
# save in TSV format
fwrite(
    agg_data[order(Sample), .(
        Sample,
        Chrom1,
        Start1,
        End1,
        Strand1,
        Chrom2,
        Start2,
        End2,
        Strand2,
        LogOdds,
        Resolution
    )],
    paste0(ARGS$prefix, ".tsv"),
    sep = "\t",
    col.names = TRUE
)

# save in BEDPE format
#   see https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format
fwrite(
    agg_data[order(Sample), .(
        Chrom1,
        Start1,
        End1,
        Chrom2,
        Start2,
        End2,
        Sample,
        LogOdds,
        Strand1,
        Strand2,
        Resolution
    )],
    paste0(ARGS$prefix, ".bedpe"),
    sep = "\t",
    col.names = FALSE
)