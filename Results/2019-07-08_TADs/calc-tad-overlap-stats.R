# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
suppressMessages(library("argparse"))

if (!interactive()) {
    PARSER <- argparse::ArgumentParser(
        description = "Calculate TAD overlap similarities between PCa samples"
    )
    ARGS <- PARSER$parse_args()
}

CHRS = paste0("chr", c(1:22, "X", "Y"))
SAMPLES = paste0("PCa", c(3023,13266,13848,14121,19121,33173,40507,51852,53687,56413,57054,57294,58215))
WINDOWS = c(30, 20, 10, 3)

# ==============================================================================
# Functions
# ==============================================================================

# ==============================================================================
# Data
# ==============================================================================
# load individual sample TAD calls
## apply over samples
domains = rbindlist(lapply(
    SAMPLES,
    function(s) {
        # apply over windows
        dt1 = rbindlist(lapply(
            WINDOWS,
            function(w) {
                ## apply over chromosomes
                dt2 = rbindlist(lapply(
                    CHRS,
                    function(chrom) {
                        dt3 = fread(
                            paste0("TADs/w_", w, "/", s, ".40000bp.", chrom, ".domains.bed"),
                            header = FALSE,
                            col.names = c("chr", "start", "end", "tag")
                        )
                        return(dt3[tag == "domain", .N, by = "chr"])
                    }
                ))
                dt2[, Window := w]
                dt2[, Resolution := 40000]
                return(dt2)
            }
        ))
        dt1[, Sample := s]
        return(dt1)
    }
))

# count TAD intersections from bedtools
## apply over pairs of samples
combos = combinations(n = length(SAMPLES), r = 2, v = SAMPLES)
intersections = rbindlist(apply(
    combos,
    1,
    function(combo) {
        # apply over windows
        print(combo)
        dt1 = rbindlist(lapply(
            WINDOWS,
            function(w) {
                ## apply over chromosomes
                fname = paste0("TADs/Comparisons/w_", w, ".", combo[1], ".", combo[2], ".intersected.bed")
                dt2 = fread(
                        fname,
                        header = FALSE,
                        col.names = c("chr", "start", "end", "tag")
                    )
                dt2[, Window := w]
                return(dt2[tag == "domain", .N, by = c("chr", "Window")])
            }
        ))
        dt1[, Sample1 := combo[1]]
        dt1[, Sample2 := combo[2]]
        return(dt1)
    }
))

# ==============================================================================
# Analysis
# ==============================================================================
# calculate percentage of similar TADs per sample
samplewise_intersections = rbindlist(lapply(
    SAMPLES,
    function(s) {
        s_ints = rbindlist(
            list(
                intersections[Sample1 == s, sum(N), by = c("Window", "Sample2")],
                intersections[Sample2 == s, sum(N), by = c("Window", "Sample1")]
            ),
            use.names = FALSE
        )
        for (w in WINDOWS) {
            sample_total = domains[Window == w & Sample == s, sum(N)]
            s_ints[Window == w, Total := sample_total]
        }
        s_ints[, Sample1 := s]
        return(s_ints)
    }
))
samplewise_intersections[, Frac := V1 / Total]
colnames(samplewise_intersections)[3] = "N_Int"
samplewise_intersections[, Window := factor(Window, levels = rev(WINDOWS), ordered = TRUE)]

# save results
fwrite(
    samplewise_intersections[, .(Sample1, Sample2, Window, N_Int, Total, Frac)],
    "TADs/Comparisons/comparison-total-counts.tsv",
    col.names = TRUE,
    sep = "\t"
)

# ==============================================================================
# Plots
# ==============================================================================
gg = (
    ggplot(data = samplewise_intersections)
    + geom_boxplot(aes(x = Window, y = 100 * Frac), outlier.shape = NA)
    + geom_point(aes(x = Window, y = 100 * Frac, colour = Sample1), position=position_jitter(width = 0.5))
    + labs(x = "Window Size", y = "Similar TADs (%)")
    + guides(colour = FALSE)
    + facet_wrap(~ Sample1)
    + theme_minimal()
    # + theme(
    #     axis.text.x = element_text(angle=60, hjust = 0.5, vjust = 1)
    # )
)
ggsave(
    "Plots/tad-similarity-counts.png",
    height = 20,
    width = 20,
    units = "cm"
)
