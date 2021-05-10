# ==============================================================================
# Meta
# ==============================================================================
# Summary breaks
# --------------------------------------
# Description: Summarizes the breakpoints and SV events in each patient
# Author: James Hawley


# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))

CHRS <- paste0("chr", c(1:22, "X", "Y"))

hg38 <- fread(
    file.path(
        "..", "..", "Data", "Processed", "2019-06-18_PCa-LowC-sequencing",
        "hg38.sizes.txt"
    ),
    sep = "\t",
    header = FALSE,
    col.names = c("Chrom", "Length")
)
hg38[, Bins := ceiling(Length / 10^6)]


# ==============================================================================
# Data
# ==============================================================================
# load metadata
metadata <- fread(
    file.path(
        "..", "..", "Data", "External", "LowC_Samples_Data_Available.tsv"
    ),
    sep = "\t",
    header = TRUE
)
metadata <- metadata[Include == "Yes"]
metadata[, SampleID := paste0("PCa", get("Sample ID"))]
SAMPLES <- metadata$SampleID

# load breakpoint data
breakpoints <- fread(
    file.path("Graphs", "sv-breakpoints.tsv"),
    sep = "\t",
    header = TRUE
)
breakpoint_pairs <- fread(
    file.path("Graphs", "sv-breakpoints.paired.tsv"),
    sep = "\t",
    header = TRUE
)


# ==============================================================================
# Analysis
# ==============================================================================
# 1. Initial data cleaning
# --------------------------------------
# convert breakpoints into Mbp bins and counts the rearrangements
breakpoints[, StartBin := floor(start / 10^6)]
breakpoints[, EndBin := floor(end / 10^6)]

# add T2E status
breakpoints <- merge(
    breakpoints,
    metadata[, .SD, .SDcols = c("SampleID", "T2E Status")],
    by = "SampleID"
)

breakpoint_pairs <- merge(
    breakpoint_pairs,
    metadata[, .SD, .SDcols = c("SampleID", "T2E Status")],
    by = "SampleID"
)


# 2. Counting breakpoints
# --------------------------------------
cat("Summarizing breakpoints\n")
# count breakpoints per sample
breakpoints_counted <- breakpoints[, .N, keyby = c("SampleID", "T2E Status")]

# convert breakpoints in genomic regions into bins
breakpoints_binned <- rbindlist(lapply(
    1:breakpoints[, .N],
    function(i) {
        start_bin <- breakpoints[i, StartBin]
        end_bin <- breakpoints[i, EndBin]
        dt <- data.table(
            SampleID = breakpoints[i, SampleID],
            chr = breakpoints[i, chr],
            Bin = start_bin:end_bin,
            Count = 1
        )
        return(dt)
    }
))
# count all breakpoints in each bin
breakpoints_summed <- breakpoints_binned[,
    sum(Count),
    by = c("SampleID", "chr", "Bin")
]
colnames(breakpoints_summed) <- c("SampleID", "chr", "Bin", "Count")
fwrite(
    breakpoints_summed[order(SampleID, chr, Bin)],
    file.path("Statistics", "breakpoints.binned.tsv"),
    sep = "\t",
    col.names = TRUE
)

# count breakpoints per chromosome
breakpoints_by_chrom <- breakpoints[,
    .N,
    keyby = c("SampleID", "T2E Status", "chr")
]
breakpoints_by_chrom[, N_per_mb := apply(.SD, 1, function(r) {
    as.numeric(r["N"]) / hg38[Chrom == r["chr"], (Length / 10^6)]
})]
fwrite(
    breakpoints_by_chrom[order(SampleID, chr)],
    file.path("Statistics", "breakpoints.by-chrom.tsv"),
    sep = "\t",
    col.names = TRUE
)


# 2. Counting breakpoint pairs
# --------------------------------------
cat("Summarizing breakpoint pairs\n")
# count the number of inter-/intra-chromosomal pairs
inter_intra_counts <- unique(breakpoint_pairs[,
    .(
        T2E_Status = get("T2E Status"),
        Total = .N,
        Intrachromosomal = sum(chr_x == chr_y),
        Interchromosomal = sum(chr_x != chr_y)
    ),
    by = "SampleID"
])
inter_intra_counts <- melt(
    inter_intra_counts,
    id.vars = c("SampleID", "T2E_Status"),
    variable.name = "Class",
    value.name = "Count"
)
fwrite(
    inter_intra_counts,
    file.path("Statistics", "breakpoint-pairs.inter-intra-chromosomal.tsv"),
    sep = "\t",
    col.names = TRUE
)


# 3. Counting SV events (simple and complex)
# --------------------------------------
cat("Summarizing SV events\n")
# calculate the length of events in each sample
# (i.e. size of each connected component)
breakpoint_components <- breakpoints[,
    .(
        N_Breakpoints = .N,
        N_Chr = length(unique(chr))
    ),
    by = c("SampleID", "T2E Status", "component_ID")
]
fwrite(
    breakpoint_components,
    file.path("Statistics", "sv-components.tsv"),
    sep = "\t",
    col.names = TRUE
)

# calculate the number of total events and complex events
breakpoint_components_counted <- breakpoint_components[,
    .(
        N_Events = .N,
        N_Complex_Events = sum(N_Breakpoints > 2)
    ),
    by = c("SampleID", "T2E Status")
]
fwrite(
    breakpoint_components_counted,
    file.path("Statistics", "sv-components.counts.tsv"),
    sep = "\t",
    col.names = TRUE
)

# 4. Hypothesis testing
# --------------------------------------
cat("Hypothesis tests\n")
# perform Mann-Whitney U to test if T2E patients have more of
# the following than non-T2E:
#   a) unique breakpoints
#   b) total events
#   c) complex events
#   d) Inter-chromosomal breakpoint pairs
#   e) Intra-chromosomal breakpoint pairs
#   f) independence between inter-/intra-chromosomal pairs and T2E status

htests <- list(
    # a)
    "breakpoints" = wilcox.test(
        x = breakpoints_counted[get("T2E Status") == "No", N],
        y = breakpoints_counted[get("T2E Status") == "Yes", N],
        alternative = "less"
    ),
    # b)
    "total" = wilcox.test(
        x = breakpoint_components_counted[get("T2E Status") == "No", N_Events],
        y = breakpoint_components_counted[get("T2E Status") == "Yes", N_Events],
        alternative = "less"
    ),
    # c)
    "complex" = wilcox.test(
        x = breakpoint_components_counted[
            get("T2E Status") == "No",
            N_Complex_Events
        ],
        y = breakpoint_components_counted[
            get("T2E Status") == "Yes",
            N_Complex_Events
        ],
        alternative = "less"
    ),
    # d)
    "inter" = wilcox.test(
        x = inter_intra_counts[
            T2E_Status == "No" & Class == "Interchromosomal",
            Count
        ],
        y = inter_intra_counts[
            T2E_Status == "Yes" & Class == "Interchromosomal",
            Count
        ],
        alternative = "less"
    ),
    # e)
    "intra" = wilcox.test(
        x = inter_intra_counts[
            T2E_Status == "No" & Class == "Intrachromosomal",
            Count
        ],
        y = inter_intra_counts[
            T2E_Status == "Yes" & Class == "Intrachromosomal",
            Count
        ],
        alternative = "less"
    ),
    # f)
    "t2e_chr" = fisher.test(
        rbind(
            c(
                inter_intra_counts[
                    Class == "Intrachromosomal" & T2E_Status == "No",
                    sum(Count)
                ],
                inter_intra_counts[
                    Class == "Intrachromosomal" & T2E_Status == "Yes",
                    sum(Count)
                ]
            ),
            c(
                inter_intra_counts[
                    Class == "Interchromosomal" & T2E_Status == "No",
                    sum(Count)
                ],
                inter_intra_counts[
                    Class == "Interchromosomal" & T2E_Status == "Yes",
                    sum(Count)
                ]
            )
        ),
        alternative = "two.sided"
    )
)
# save hypothesis test results
saveRDS(htests, file.path("Statistics", "htests.rds"))

cat("----\n\tTotal breakpoints:\n")
print(htests$breakpoints)
cat("----\n\tTotal events:\n")
print(htests$total)
cat("----\n\tComplex events:\n")
print(htests$complex)
cat("----\n\tInter-chromosomal pairs:\n")
print(htests$inter)
cat("----\n\tIntra-chromosomal pairs:\n")
print(htests$intra)
cat("----\n\tBreakpoint pair location and T2E status:\n")
print(htests$t2e_chr)