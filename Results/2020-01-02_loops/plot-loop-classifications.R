# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))

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
# count loop calls per overlapping CRE type
loop_cre_counts = loops[, .N, by = c("SampleID", "CRE_Overlap")]
loop_cre_counts[, Frequency := 0]

# calculate frequencies, not just total counts
invisible(sapply(
    1:loop_cre_counts[, .N],
    function(i) {
        count = loop_cre_counts[i, N]
        total = n_loops[SampleID == loop_cre_counts[i, SampleID], N]
        loop_cre_counts[i, Frequency := count / total]
}))

# count loop calls per overlapping CRE type
loop_bound_counts = loops[CRE_Overlap == 0, .N, by = c("SampleID", "TAD_Overlap")]
loop_bound_counts[, Frequency := 0]

# calculate frequencies, not just total counts
invisible(sapply(
    1:loop_bound_counts[, .N],
    function(i) {
        count = loop_bound_counts[i, N]
        total = n_loops[SampleID == loop_bound_counts[i, SampleID], N]
        loop_bound_counts[i, Frequency := count / total]
}))

# ==============================================================================
# Plots
# ==============================================================================
gg = (
    ggplot(data = loop_cre_counts)
    + geom_col(
        aes(x = SampleID, y = N, fill = factor(CRE_Overlap, levels = 0:2, ordered = TRUE)),
        position = "stack"
    )
    + labs(x = NULL, y = "Number of Loops")
    + scale_fill_manual(
        limits = 0:2,
        breaks = 0:2,
        values = c(
            "#ece7f2",
            "#a6bddb",
            "#2b8cbe"
        )
    )
    + theme_minimal()
    + theme(
        axis.text.x = element_text(angle = 90)
    )
)
ggsave(
    file.path("Plots", "loop-CRE-overlap.count.png"),
    height = 12,
    width = 20,
    units = "cm"
)

gg = (
    ggplot(data = loop_cre_counts)
    + geom_col(
        aes(x = SampleID, y = N, fill = factor(CRE_Overlap, levels = 0:2, ordered = TRUE)),
        position = "fill"
    )
    + labs(x = NULL, y = "Frequency of Loops")
    + scale_fill_manual(
        limits = 0:2,
        breaks = 0:2,
        values = c(
            "#ece7f2",
            "#a6bddb",
            "#2b8cbe"
        )
    )
    + theme_minimal()
    + theme(
        axis.text.x = element_text(angle = 90)
    )
)
ggsave(
    file.path("Plots", "loop-CRE-overlap.proportion.png"),
    height = 12,
    width = 20,
    units = "cm"
)

gg = (
    ggplot(data = loop_bound_counts)
    + geom_col(
        aes(x = SampleID, y = N, fill = factor(TAD_Overlap, levels = 0:2, ordered = TRUE)),
        position = "stack"
    )
    + labs(x = NULL, y = "Number of Loops")
    + scale_fill_manual(
        limits = 0:2,
        breaks = 0:2,
        values = c(
            "#fff7bc",
            "#fec44f",
            "#d95f0e"
        )
    )
    + theme_minimal()
    + theme(
        axis.text.x = element_text(angle = 90)
    )
)
ggsave(
    file.path("Plots", "loop-TAD-overlap.count.png"),
    height = 12,
    width = 20,
    units = "cm"
)

gg = (
    ggplot(data = loop_bound_counts)
    + geom_col(
        aes(x = SampleID, y = N, fill = factor(TAD_Overlap, levels = 0:2, ordered = TRUE)),
        position = "fill"
    )
    + labs(x = NULL, y = "Frequency of Loops")
    + scale_fill_manual(
        limits = 0:2,
        breaks = 0:2,
        values = c(
            "#fff7bc",
            "#fec44f",
            "#d95f0e"
        )
    )
    + theme_minimal()
    + theme(
        axis.text.x = element_text(angle = 90)
    )
)
ggsave(
    file.path("Plots", "loop-TAD-overlap.proportion.png"),
    height = 12,
    width = 20,
    units = "cm"
)

# ==============================================================================
# Save data
# ==============================================================================
fwrite(
    loop_cre_counts,
    file.path("Classification", "loops.CREs.tsv"),
    sep = "\t",
    col.names = TRUE
)

fwrite(
    loop_bound_counts,
    file.path("Classification", "loops.TADs.tsv"),
    sep = "\t",
    col.names = TRUE
)
