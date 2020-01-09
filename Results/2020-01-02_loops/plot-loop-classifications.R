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
# load sample metadata
metadata = fread(file.path("..", "..", "Data", "External", "LowC_Samples_Data_Available.tsv"))
colnames(metadata) = gsub(" ", "_", colnames(metadata))
metadata[, Sample_ID := paste0("PCa", Sample_ID)]

# keep track of the number of loops called per patient
n_loops = data.table(
    SampleID = metadata[, Sample_ID],
    N = 0
)

loops = rbindlist(lapply(
    metadata[, Sample_ID],
    function(id) {
        dt = fread(file.path("Classification", paste0(id, ".loops.tsv")))
        dt[, SampleID := id]
        n_loops[SampleID == id, N := dt[, .N]]
        return(dt)
    }
))


# ==============================================================================
# Analysis
# ==============================================================================
# count loop calls per overlapping type
loop_counts = loops[, .N, by = c("SampleID", "CRE_Overlap")]
loop_counts[, Frequency := 0]

# calculate frequencies, not just total counts
invisible(sapply(
    1:loop_counts[, .N],
    function(i) {
        count = loop_counts[i, N]
        total = n_loops[SampleID == loop_counts[i, SampleID], N]
        loop_counts[i, Frequency := count / total]
}))

# combine `notx_and_y` and `x_and_noty` together, so the totals are number of loops overlapping 0, 1, or 2 CREs
loop_counts_by_num_peaks = data.table(
    SampleID = rep(metadata[, Sample_ID], 3),
    Overlapping_CREs = rep(0:2, each = metadata[, .N]),
    N = 0,
    Frequency = 0
)
# calculate new frequencies and counts
invisible(sapply(
    1:loop_counts_by_num_peaks[, .N],
    function(i) {
        s = loop_counts_by_num_peaks[i, SampleID]
        n_overlaps = loop_counts_by_num_peaks[i, Overlapping_CREs]
        if (n_overlaps == 0) {
            count = loop_counts[SampleID == s & CRE_Overlap == "notx_and_noty", N]
        } else if (n_overlaps == 1) {
            count = loop_counts[SampleID == s & CRE_Overlap %in% c("notx_and_y", "x_and_noty"), sum(N)]
        } else {
            count = loop_counts[SampleID == s & CRE_Overlap == "x_and_y", N]
        }
        total = n_loops[SampleID == loop_counts_by_num_peaks[i, SampleID], N]
        loop_counts_by_num_peaks[i, N := count]
        loop_counts_by_num_peaks[i, Frequency := count / total]
}))


# ==============================================================================
# Plots
# ==============================================================================
gg = (
    ggplot(data = loop_counts)
    + geom_col(aes(x = SampleID, y = N, fill = CRE_Overlap), position = "stack")
    + labs(x = NULL, y = "Number of Loops")
    + theme_minimal()
    + theme(
        axis.text.x = element_text(angle = 90)
    )
)
ggsave(
    file.path("Plots", "loop-CRE-overlap.all.count.png"),
    height = 12,
    width = 20,
    units = "cm"
)

gg = (
    ggplot(data = loop_counts)
    + geom_col(aes(x = SampleID, y = N, fill = CRE_Overlap), position = "fill")
    + labs(x = NULL, y = "Frequency of Loops")
    + theme_minimal()
    + theme(
        axis.text.x = element_text(angle = 90)
    )
)
ggsave(
    file.path("Plots", "loop-CRE-overlap.all.proportion.png"),
    height = 12,
    width = 20,
    units = "cm"
)

gg = (
    ggplot(data = loop_counts_by_num_peaks)
    + geom_col(aes(x = SampleID, y = N, fill = Overlapping_CREs), position = "stack")
    + labs(x = NULL, y = "Number of Loops")
    + theme_minimal()
    + theme(
        axis.text.x = element_text(angle = 90)
    )
)
ggsave(
    file.path("Plots", "loop-CRE-overlap.sum.count.png"),
    height = 12,
    width = 20,
    units = "cm"
)

gg = (
    ggplot(data = loop_counts_by_num_peaks)
    + geom_col(aes(x = SampleID, y = N, fill = Overlapping_CREs), position = "fill")
    + labs(x = NULL, y = "Frequency of Loops")
    + theme_minimal()
    + theme(
        axis.text.x = element_text(angle = 90)
    )
)
ggsave(
    file.path("Plots", "loop-CRE-overlap.sum.proportion.png"),
    height = 12,
    width = 20,
    units = "cm"
)

# ==============================================================================
# Save data
# ==============================================================================
fwrite(
    loop_counts,
    file.path("Classification", "all.loops.all.tsv"),
    sep = "\t",
    col.names = TRUE
)

fwrite(
    loop_counts_by_num_peaks,
    file.path("Classification", "all.loops.sum.tsv"),
    sep = "\t",
    col.names = TRUE
)
