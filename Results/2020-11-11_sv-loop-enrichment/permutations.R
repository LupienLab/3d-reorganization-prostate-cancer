# ==============================================================================
# Meta
# ==============================================================================
# permutations
# --------------------------------------
# Description: Test for enrichment of loops in SV breakpoints based on whether they affect expression or not
# Author: James Hawley


# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("logging"))

loginfo("Loading packages")
suppressMessages(library("data.table"))
suppressMessages(library("GenomicRanges"))

# ==============================================================================
# Functions
# ==============================================================================

# ==============================================================================
# Data
# ==============================================================================
loginfo("Loading data")

# load metadata
metadata <- fread("config.tsv", sep = "\t")

# load SV breakpoints and their expression changes
sv_bp <- fread("../2020-02-19_chromoplexy/Graphs/sv-breakpoints.tsv", sep = "\t")
sv_exprs <- fread("../2020-06-18_sv-disruption-expression/1sample-results.tsv", sep = "\t")

# load loop anchors
loops <- fread("../2020-10-06_loops/Loops/merged-loops.sample-counts.tsv", sep = "\t")
anchors <- rbind(
    loops[, .(
        chr = chr_x,
        start = start_x,
        end = end_x,
        anchor_ID = anchor_ID_x,
        loop_ID = loop_ID,
        fdr = fdr,
        detection_scale = detection_scale
    )],
    loops[, .(
        chr = chr_y,
        start = start_y,
        end = end_y,
        anchor_ID = anchor_ID_y,
        loop_ID = loop_ID,
        fdr = fdr,
        detection_scale = detection_scale
    )]
)

# ==============================================================================
# Analysis
# ==============================================================================
loginfo("Calculating observed intersection of SV breakpoints and loops")

# merge SV breakpoints and expression changes
any_dge <- sv_exprs[Significant == TRUE, .(DGE = sum(N) > 0), by = "test_ID"]
sv_bp <- merge(
    x = sv_bp,
    y = any_dge,
    by = "test_ID"
)

# convert anchors and SVs to GRanges
anch_gr <- GRanges(
    seqnames = anchors$chr,
    ranges = IRanges(
        start = anchors$start + 1,
        end = anchors$end
    ),
    anchor_ID = anchors$anchor_ID,
    loop_ID = anchors$loop_ID,
    fdr = anchors$fdr,
    detection_scale = anchors$detection_scale
)

svbp_gr <- GRanges(
    seqnames = sv_bp$chr,
    ranges = IRanges(
        start = sv_bp$start,
        end = sv_bp$end
    ),
    breakpoint_ID = sv_bp$breakpoint_ID,
    component_ID = sv_bp$component_ID,
    test_ID = sv_bp$test_ID
)

# find overlaps between anchors and SV breakpoints
hits <- as.data.table(findOverlaps(anch_gr, svbp_gr))

# calculate number of overlaps by SV (i.e. component_ID)
sv_n_loops <- cbind(1, sv_bp[hits$subjectHits])
sv_n_loops <- sv_n_loops[,
    .(
        n_loops = sum(V1),
        any_dge = any(DGE)
    ),
    keyby = c("SampleID", "component_ID")
]

obs_enrichment <- sv_n_loops[, .(Median_Loops = median(n_loops)), keyby = DGE]
obs_log2fc <- (
    obs_enrichment[DGE == TRUE, log2(Median_Loops)]
    - obs_enrichment[DGE == FALSE, log2(Median_Loops)]
)


loginfo("Calculating permutations")
n_dge <- sv_n_loops[any_dge == TRUE, .N]
n_total <- sv_n_loops[, .N]

set.seed(1111)
n_perms <- 100000
perm_fcs <- data.table(
    Iteration = 1:n_perms,
    median_loops_dge = 0,
    median_loops_nondge = 0
)
for (i in 1:n_perms) {
    if (i %% 1000 == 0) {
        cat(i, "\n")
    }
    perm_intersect <- copy(sv_n_loops)
    # pick n_dge SVs to count as altering expression
    idx <- sample.int(n = n_total, size = n_dge)
    perm_intersect[idx, any_dge := TRUE]
    perm_intersect[-idx, any_dge := FALSE]
    perm_enrichment <- perm_intersect[,
        .(Median_Loops = median(n_loops)),
        keyby = any_dge
    ]
    n_loops_dge <- perm_enrichment[any_dge == TRUE, Median_Loops]
    n_loops_nondge <- perm_enrichment[any_dge == FALSE, Median_Loops]
    perm_fcs[i, median_loops_dge := n_loops_dge]
    perm_fcs[i, median_loops_nondge := n_loops_nondge]
}

# calculate fold change for all permutations
perm_fcs[, dge_log2Fold := log2(median_loops_dge / median_loops_nondge)]

# calculate p-value from permutations
# number of permutations with a median fold change in number of loops based on SVs affecting expression
pval <- perm_fcs[dge_log2Fold >= obs_log2fc, .N] / perm_fcs[, .N]
cat("Permutation p-value: ", pval, "\n")

# ==============================================================================
# Save data
# ==============================================================================
loginfo("Saving data")

fwrite(
    sv_n_loops,
    "sv-loop-intersection.tsv",
    sep = "\t",
    col.names = TRUE
)

fwrite(
    obs_enrichment,
    "sv-loop-intersection.observed-enrichment.tsv",
    sep = "\t",
    col.names = TRUE
)

fwrite(
    perm_fcs,
    "sv-loop-intersection.permutations.tsv",
    sep = "\t",
    col.names = TRUE
)