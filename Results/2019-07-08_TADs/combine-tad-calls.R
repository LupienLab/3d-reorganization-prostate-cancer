# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))

# ==============================================================================
# Data
# ==============================================================================
SAMPLES = paste0('PCa', c(
    3023,
    13266,
    13848,
    14121,
    19121,
    33173,
    40507,
    51852,
    53687,
    56413,
    57054,
    57294,
    58215
))

CHRS = paste0("chr", c(1:22, "X", "Y"))

tad_stats = data.table(
    Sample = rep(SAMPLES, each = length(3:30)),
    Resolution = 40000,
    w = rep(3:30, 13),
    N = 0
)

# ==============================================================================
# Analysis
# ==============================================================================
for (i in 1:tad_stats[, .N]) {
    cat(i, "\n")
    sample_name = tad_stats[i, Sample]
    res = tad_stats[i, Resolution]
    w = tad_stats[i, w]
    tad_calls = rbindlist(lapply(
        CHRS,
        function(chrom) {
            dt = fread(
                paste0("TADs/w_", w, "/", sample_name, ".", res, "bp.", chrom, ".domains.tsv"),
                sep = "\t",
                header = TRUE
            )
            dt[, Sample := sample_name]
            dt[, w := w]
            dt[, Resolution := res]
            return(dt)
        }
    ))
    # append data to output file
    fwrite(
        tad_calls,
        "TADs/tad-calls.tsv",
        sep = "\t",
        col.names = TRUE,
        append = TRUE
    )
    tad_stats[i, N := tad_calls[tag == "domain", .N]]
    tad_stats[i, Mean_Size := tad_calls[tag == "domain", mean(to.coord - from.coord)]]
    tad_stats[i, SD_Size := tad_calls[tag == "domain", sd(to.coord - from.coord)]]
}

# save statistics on TAD calls
fwrite(
    tad_stats,
    "TADs/tad-call-stats.tsv",
    sep = "\t",
    col.names = TRUE
)
