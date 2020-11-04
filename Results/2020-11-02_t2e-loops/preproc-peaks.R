# ==============================================================================
# Meta
# ==============================================================================
# preprocess-loops
# --------------------------------------
# Description: Produce a subset of loops that overlap gene promoters
# Author: James Hawley


# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("logging"))

loginfo("Loading packages")
suppressMessages(library("data.table"))
suppressMessages(library("regioneR"))

# ==============================================================================
# Functions
# ==============================================================================

# ==============================================================================
# Data
# ==============================================================================
loginfo("Loading data")
promoters <- fread(
    "../../Data/External/GENCODE/gencode.v33.all-genes.promoters.bed",
    sep = "\t",
    header = FALSE,
    col.names = c("chr", "start", "end", "strand", "gene_id", "gene_name")
)

loops <- fread(
    "loops.overlapping-promoters.tsv",
    sep = "\t",
    header = TRUE
)

anchors <- rbind(
    loops[, .(
        chr = chr_x,
        start = start_x,
        end = end_x,
        anchor_ID = anchor_ID_x,
        loop_ID = loopID
    )],
    loops[, .(
        chr = chr_y,
        start = start_y,
        end = end_y,
        anchor_ID = anchor_ID_y,
        loop_ID = loopID
    )]
)

peaks <- fread(
    file.path("..", "2020-06-12_sv-disruption-acetylation", "Acetylation", "T2E", "t2e.all.tsv"),
    sep = "\t",
    header = TRUE
)
peaks[, peak_ID := .I]


# ==============================================================================
# Analysis
# ==============================================================================
loginfo("Intersecting regions")

prom_gr <- toGRanges(promoters)
anch_gr <- toGRanges(anchors)
peak_gr <- toGRanges(peaks)

# exclude peaks that overlap promoters
prom_hits <- findOverlaps(peak_gr, prom_gr)
peak_gr <- peak_gr[-queryHits(prom_hits)]

# only include peaks that overlap an anchor
anch_hits <- findOverlaps(peak_gr, anch_gr)
overlapping_peaks <- peak_gr[queryHits(anch_hits)]
ov_peaks_dt <- peaks[peak_ID %in% overlapping_peaks$peak_ID]

# only include loops that overlap at least one enhancer, too
overlapping_anchs <- anch_gr[subjectHits(anch_hits)]
ov_loops_dt <- loops[loopID %in% overlapping_anchs$loop_ID]

# ==============================================================================
# Save data
# ==============================================================================
loginfo("Saving data")

fwrite(
    ov_loops_dt,
    "loops.overlapping-promoters.overlapping-enhancers.tsv",
    sep = "\t",
    col.names = TRUE
)
fwrite(
    ov_peaks_dt[,
        .SD,
        keyby = c("chr", "start", "end"),
        .SDcols = c("chr", "start", "end", "Conc", "Conc_T2E", "Conc_NonT2E", "Fold", "FDR", "peak_ID")
    ],
    "enhancers.overlapping-loops.tsv",
    sep = "\t",
    col.names = TRUE
)
