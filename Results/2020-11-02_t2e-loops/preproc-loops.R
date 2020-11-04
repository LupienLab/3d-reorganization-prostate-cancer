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
    file.path("..", "2020-10-06_loops", "Loops", "merged-loops.sample-counts.tsv"),
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

# ==============================================================================
# Analysis
# ==============================================================================
loginfo("Intersecting regions")

prom_gr <- toGRanges(promoters)
anch_gr <- toGRanges(anchors)

hits <- findOverlaps(anch_gr, prom_gr)
overlapping_anchors <- anch_gr[queryHits(hits)]
ov_loops_dt <- loops[loopID %in% overlapping_anchors$loop_ID]

overlapping_genes <- prom_gr[subjectHits(hits)]
ov_genes_dt <- promoters[gene_id %in% overlapping_genes$gene_id]

# ==============================================================================
# Save data
# ==============================================================================
loginfo("Saving data")

fwrite(
    ov_loops_dt,
    "loops.overlapping-promoters.tsv",
    sep = "\t",
    col.names = TRUE
)
fwrite(
    ov_genes_dt,
    "promoters.overlapping-loops.tsv",
    sep = "\t",
    col.names = TRUE
)
