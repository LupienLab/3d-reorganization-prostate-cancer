# ==============================================================================
# Meta
# ==============================================================================
# acetyl-loop-exprs
# --------------------------------------
# Description: Determine the relationship between differential acetylation, loops, and expression in T2E+/- samples
# Author: James Hawley


# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("logging"))

loginfo("Loading packages")
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
suppressMessages(library("regioneR"))

# ==============================================================================
# Functions
# ==============================================================================

# ==============================================================================
# Data
# ==============================================================================
loginfo("Loading data")
intersected_loop_acetyl <- fread("loops.intersected-sig-acetyl.tsv", sep = "\t")

# differential expression results
so_transcripts <- fread("../2020-06-18_sv-disruption-expression/results.transcripts.tsv")
so_transcripts <- so_transcripts[complete.cases(so_transcripts)]
so_genes <- fread("../2020-06-18_sv-disruption-expression/results.genes.tsv")

# GENCODE annotations
gencode_genes <- fread("../../Data/External/GENCODE/gencode.v33.all-genes.bed")
gencode_tx <- fread("../../Data/External/GENCODE/gencode.v33.all-transcripts.bed")

# ==============================================================================
# Analysis
# ==============================================================================
loginfo("Performing calculations")
important_loops <- intersected_loop_acetyl[peak_FDR < 0.05 & abs(peak_Fold) > 1 & loop_Loop_Type != "shared"]
important_anchors <- unique(rbind(
    important_loops[, .(
        chr = loop_chr_x,
        start = loop_start_x,
        end = loop_end_x,
        anchor_ID = loop_anchor_ID_x,
        loop_ID = loop_loopID,
        loop_type = ifelse(loop_Loop_Type == "T2E-specific", "T2E+", "T2E-"),
        peak_type = ifelse(peak_Fold > 1, "T2E+", "T2E-")
    )],
    important_loops[, .(
        chr = loop_chr_y,
        start = loop_start_y,
        end = loop_end_y,
        anchor_ID = loop_anchor_ID_y,
        loop_ID = loop_loopID,
        loop_type = ifelse(loop_Loop_Type == "T2E-specific", "T2E+", "T2E-"),
        peak_type = ifelse(peak_Fold > 1, "T2E+", "T2E-")
    )]
))

# remove loops where the is both a peak gained and lost in an anchor in the T2E+ samples
# this is to remove ambiguity of the effect of acetylation
important_anchors[, .N, by = c("loop_ID", "anchor_ID")][N > 1, length(unique(loop_ID))] # there are only 2 occurrences of this
ambiguous_loops_to_remove <- important_anchors[, .N, by = c("loop_ID", "anchor_ID")][N > 1, unique(loop_ID)]
# make a "not in" convenience function
"%ni%" <- Negate("%in%")
# remove these ambiguous loops
important_anchors <- important_anchors[loop_ID %ni% ambiguous_loops_to_remove]


important_regions <- unique(c(
    toGRanges(
        important_loops[, .(
            chr = loop_chr_x,
            start = loop_start_x,
            end = loop_end_x,
            loop_ID = loop_loopID,
            loop_type = ifelse(loop_Loop_Type == "T2E-specific", "T2E+", "T2E-"),
            peak_type = ifelse(peak_Fold > 1, "T2E+", "T2E-")
        )]
    ),
    toGRanges(
        important_loops[, .(
            chr = loop_chr_y,
            start = loop_start_y,
            end = loop_end_y,
            loop_ID = loop_loopID,
            loop_type = ifelse(loop_Loop_Type == "T2E-specific", "T2E+", "T2E-"),
            peak_type = ifelse(peak_Fold > 1, "T2E+", "T2E-")
        )]
    )
))

# find all transcripts that overlap the important regions (i.e. anchors of any loop that also contains differential acetylation)
tx_gr <- toGRanges(so_transcripts)
hits <- as.data.table(findOverlaps(tx_gr, important_regions))

# map hits back to the transcripts
# remove and long transcripts that overlap with multiple anchors so that
# each transcript, if it overlaps a loop, will only overlap a single region
hits[, .N, by = "queryHits"][N > 1, .N] # only 1 case for PEBP4-201 transcript
# remove the transcript
query_hit_to_remove <- hits[, .N, by = "queryHits"][N > 1, queryHits]
hits <- hits[queryHits != query_hit_to_remove]

# merge in transcript information to overlaps
anchor_exprs <- merge(
    x = so_transcripts[, .(.SD, .I)],
    y = hits,
    by.x = "I",
    by.y = "queryHits",
    all.x = FALSE,
    all.y = TRUE
)
colnames(anchor_exprs) <- gsub("^.SD.", "", colnames(anchor_exprs))

# merge in loop information to overlaps
anchor_exprs <- merge(
    x = anchor_exprs,
    y = 
)

# ==============================================================================
# Save data
# ==============================================================================
loginfo("Saving data")