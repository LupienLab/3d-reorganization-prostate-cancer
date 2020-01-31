# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
suppressMessages(library("gtools"))


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
                    paste0("../2020-01-15_TAD-aggregation/resolved-TADs/separated-TADs/", s, ".40000bp.w_", w, ".domains.bed"),
                        header = FALSE,
                        col.names = c("chr", "start", "end", "lower_persistence", "upper_persistence")
                    )
                return(dt2[, .(w, .N), by = "chr"])
            }
        ))
        colnames(dt1) = c("chr", "w", "N")
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
                fname = paste0("TAD-comparisons/intersected/w_", w, ".", combo[1], ".", combo[2], ".intersected.bed")
                dt2 = fread(
                        fname,
                        header = FALSE,
                        col.names = c("chr", "start", "end", "lower_persistence", "upper_persistence")
                    )
                dt2[, Window := w]
                return(dt2[, .N, by = c("chr", "Window")])
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
                # all combos where its sample 1
                intersections[Sample1 == s, sum(N), by = c("Window", "Sample2")],
                # all combos where its sample 2
                intersections[Sample2 == s, sum(N), by = c("Window", "Sample1")]
            ),
            use.names = FALSE
        )
        # calculate the number of domains for that sample
        for (w_size in WINDOWS) {
            sample_total = domains[w == w_size & Sample == s, sum(N)]
            s_ints[Window == w_size, Total1 := sample_total]
        }
        # add sample name as a column
        s_ints[, Sample1 := s]
        return(s_ints[, .SD, .SDcols = c("Sample1", "Sample2", "Window", "V1", "Total1")])
    }
))
colnames(samplewise_intersections)[4] = "N_Int"
samplewise_intersections[, Frac := N_Int / Total1]

# ==============================================================================
# Save data
# ==============================================================================
fwrite(
    samplewise_intersections,
    "TAD-comparisons/comparison-total-counts.tsv",
    col.names = TRUE,
    sep = "\t"
)
