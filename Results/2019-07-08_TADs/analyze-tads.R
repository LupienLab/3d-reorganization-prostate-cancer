# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))


# ==============================================================================
# Functions
# ==============================================================================

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
    Sample = rep(SAMPLES, each = 3),
    Resolution = 40000,
    w = rep(c(10, 20, 30), 13),
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
            return(dt)
        }
    ))
    tad_stats[i, N := tad_calls[tag == "domain", .N]]
    tad_stats[i, Mean_Size := tad_calls[tag == "domain", mean(to.coord - from.coord)]]
    tad_stats[i, SD_Size := tad_calls[tag == "domain", sd(to.coord - from.coord)]]
}

fwrite(
    tad_stats,
    "TADs/tad-call-stats.tsv",
    sep = "\t",
    col.names = TRUE
)

# ==============================================================================
# Plots
# ==============================================================================
gg = (
    ggplot(data = tad_stats)
    + geom_point(aes(x = w, y = N))
    + labs(x = "Window Size", y = "Number of TADs")
    + facet_wrap(~ Sample)
)
ggsave(
    "Plots/tad-counts.png",
    height = 12,
    width = 20,
    units = "cm"
)

gg = (
    ggplot(data = tad_stats)
    + geom_pointrange(aes(
        x = w,
        y = Mean_Size,
        ymin = Mean_Size - SD_Size / sqrt(N),
        ymax = Mean_Size + SD_Size / sqrt(N)
    ))
    + labs(x = "Window Size", y = "Size of TADs")
    + facet_wrap(~ Sample)
)
ggsave(
    "Plots/tad-sizes.png",
    height = 12,
    width = 20,
    units = "cm"
)
