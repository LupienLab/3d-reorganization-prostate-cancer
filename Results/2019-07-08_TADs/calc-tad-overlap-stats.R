# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
suppressMessages(library("argparse"))
suppressMessages(library("gtools"))

if (!interactive()) {
    PARSER <- argparse::ArgumentParser(
        description = "Calculate TAD overlap similarities between PCa samples"
    )
    ARGS <- PARSER$parse_args()
}

# ==============================================================================
# Data
# ==============================================================================
# load metadata
metadata = fread("../../Data/External/LowC_Samples_Data_Available.tsv")
colnames(metadata) = gsub(" ", "_", colnames(metadata))

SAMPLES = metadata[, paste0("PCa", Sample_ID)]
CHRS = paste0("chr", c(1:22, "X", "Y"))
WINDOWS = 3:30

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
                dt2 = fread(
                    paste0("TADs/w_", w, "/", s, ".40000bp.domains.bed"),
                        header = FALSE,
                        col.names = c("chr", "start", "end", "tag")
                    )
                return(dt2[tag == "domain", .(w, 40000, .N), by = "chr"])
            }
        ))
        colnames(dt1) = c("chr", "w", "res", "N")
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
        for (w_size in WINDOWS) {
            sample_total = domains[w == w_size & Sample == s, sum(N)]
            s_ints[Window == w_size, Total := sample_total]
        }
        s_ints[, Sample1 := s]
        return(s_ints)
    }
))
samplewise_intersections[, Frac := V1 / Total]
colnames(samplewise_intersections)[3] = "N_Int"
samplewise_intersections[, Window := factor(Window, levels = WINDOWS, ordered = TRUE)]

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
    # + geom_point(aes(x = Window, y = 100 * Frac, colour = Sample1), position=position_jitter(width = 0.5))
    + labs(x = "Window Size", y = "Similar TADs (%)")
    + guides(colour = FALSE)
    + facet_wrap(~ Sample1)
    + theme_minimal()
    + theme(
        axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 1)
    )
)
ggsave(
    "Plots/tad-similarity-counts.png",
    height = 20,
    width = 20,
    units = "cm"
)
