# ==============================================================================
# Meta
# ==============================================================================
# preprocess-loops
# --------------------------------------
# Description: Produce a subset of loops that overlap genes
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
#' Reduce the number of genes, loops, and enhancers being considered
#'
#' @param gnbds GRanges of gene bodies
#' @param anchs GRanges of loop anchors
#' @param enhns GRanges of enhancers
#' @return returns
reduce_elements <- function(gnbds, anchs, enhns) {
    # only consider loops where an anchor overlaps a gene body and an enhancer
    anch_gnbd_hits <- findOverlaps(anchs, gnbds)
    anch_enhn_hits <- findOverlaps(anchs, enhns)
    reduced_loop_IDs <- union(
        anchs[queryHits(anch_gnbd_hits)]$loop_ID,
        anchs[queryHits(anch_enhn_hits)]$loop_ID
    )
    reduced_anch <- anchs[anchs$loop_ID %in% reduced_loop_IDs]
    
    # only consider genes that overlap a loop anchor
    reduced_gnbd <- subsetByOverlaps(gnbds, reduced_anch)
    
    # only consider enhancers that overlap a loop anchor
    reduced_enhn <- subsetByOverlaps(enhns, reduced_anch)
    return(list(
        "gnbd" = reduced_gnbd,
        "enhn" = reduced_enhn,
        "anch" = reduced_anch
    ))
}

# ==============================================================================
# Data
# ==============================================================================
loginfo("Loading data")
genes <- fread(
    "../../Data/External/GENCODE/gencode.v33.all-genes.bed",
    sep = "\t",
    header = FALSE,
    col.names = c("chr", "start", "end", "strand", "gene_id", "gene_name")
)

loops <- rbindlist(lapply(
    c("nonT2E-specific", "T2E-specific", "shared"),
    function(loop_type) {
        dt <- fread(
            paste0("loops.", loop_type, ".tsv"),
            sep = "\t",
            header = TRUE
        )
        dt[, loop_type := loop_type]
        return(dt)
    }
))

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

# coerce to GRanges for intersection calculations
gene_gr <- toGRanges(genes)
anch_gr <- toGRanges(anchors)
peak_gr <- toGRanges(peaks)

# exclude H3K27ac peaks that overlap genes
enh_prom_hits <- findOverlaps(peak_gr, gene_gr)
enhn_gr <- peak_gr[-queryHits(enh_prom_hits)]

# repeat reduction process until no changes to loops, genes, and enhancers occurs
counter <- 0
while (TRUE) {
    counter <- counter + 1
    print(counter)
    reduction <- reduce_elements(gene_gr, anch_gr, enhn_gr)
    reduced_gnbd <- reduction$gnbd
    reduced_anch <- reduction$anch
    reduced_enhn <- reduction$enhn
    if (all(identical(reduced_gnbd, gene_gr), identical(reduced_anch, anch_gr), identical(reduced_enhn, enhn_gr))) {
        break
    } else {
        gene_gr <- reduced_gnbd
        anch_gr <- reduced_anch
        enhn_gr <- reduced_enhn
    }
}

ov_loops_dt <- loops[loopID %in% reduced_anch$loop_ID]
ov_genes_dt <- genes[gene_id %in% reduced_gnbd$gene_id]
ov_enhns_dt <- peaks[peak_ID %in% reduced_enhn$peak_ID]

# ==============================================================================
# Save data
# ==============================================================================
loginfo("Saving data")

fwrite(
    ov_loops_dt,
    "overlapping.loops.tsv",
    sep = "\t",
    col.names = TRUE
)
fwrite(
    ov_genes_dt,
    "overlapping.genes.tsv",
    sep = "\t",
    col.names = TRUE
)
fwrite(
    ov_enhns_dt,
    "overlapping.enhancers.tsv",
    sep = "\t",
    col.names = TRUE
)
