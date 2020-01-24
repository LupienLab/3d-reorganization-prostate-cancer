# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("rtracklayer"))
suppressMessages(library("GenomicRanges"))
suppressMessages(library("argparse"))

if (!interactive()) {
    PARSER <- argparse::ArgumentParser(
        description = "Liftover a set of coordinates from hg19 to hg38"
    )
    PARSER$add_argument(
        "input",
        type = "character",
        help = "Input BED file"
    )
    PARSER$add_argument(
        "output",
        type = "character",
        help = "Output BED file"
    )
    ARGS <- PARSER$parse_args()
}

# ==============================================================================
# Functions
# ==============================================================================

# ==============================================================================
# Data
# ==============================================================================
# read input file
input = fread(ARGS$input, header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "ID"))
# add 1 to convert from 0-indexed positions
input[, start := start + 1]

# convert to GRanges
gr = GRanges(input)

# read liftover chain
chain = import.chain("../../hg19ToHg38.over.chain")

# ==============================================================================
# Analysis
# ==============================================================================
# lift over to hg38
gr_hg38 = liftOver(gr, chain)

# lifted over regions are separated by a few bases, so merge them together
gr_hg38_resolved = unlist(GRangesList(lapply(
    gr_hg38,
    function(r) {
        new_gr = range(r)
        new_gr$ID = r[1]$ID
        return(new_gr)
    }
)))

# ==============================================================================
# Save data
# ==============================================================================
gr_hg38_resolved_dt = data.table(
    chr = as.character(seqnames(gr_hg38_resolved)),
    start = start(gr_hg38_resolved) - 1,
    end = start(gr_hg38_resolved),
    ID = gr_hg38_resolved$ID
)

fwrite(gr_hg38_resolved_dt, ARGS$output, sep = "\t", col.names = FALSE)
