# ==============================================================================
# Meta
# ==============================================================================
# diff-acetyl-diff-loop
# --------------------------------------
# Description: Analysis of the T2E+/- specific loops and acetylation
# Author: James Hawley


# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("logging"))

loginfo("Loading packages")
suppressMessages(library("data.table"))
suppressMessages(library("GenomicRanges"))
suppressMessages(library("regioneR"))
suppressMessages(library("ggplot2"))


LOOP_TYPES <- c("shared", "T2E-specific", "nonT2E-specific")


# ==============================================================================
# Functions
# ==============================================================================

# ==============================================================================
# Data
# ==============================================================================
loginfo("Loading data")
# load metadata
metadata <- fread("config.tsv", sep = "\t")
metadata <- metadata[Include == "Yes"]
SAMPLES <- metadata[, SampleID]

# load tumour loops
loops <- rbindlist(lapply(
    LOOP_TYPES,
    function(t) {
        dt <- fread(
            paste0("loops.", t, ".tsv"),
            sep = "\t",
            header = TRUE
        )
        dt[, Loop_Type := t]
        return(dt)
    }
))

# load T2E+/- differential acetylation results
acetyl <- fread(
    "../2020-06-12_sv-disruption-acetylation/Acetylation/T2E/t2e.all.tsv",
    sep = "\t",
    header = TRUE
)

acetyl_sig <- acetyl[FDR < 0.05]

# ==============================================================================
# Analysis
# ==============================================================================
loginfo("Overlapping loops with peaks")

# convert loops to anchors for overlap
anchors <- unique(rbind(
    loops[, .(
        chr = chr_x,
        start = start_x,
        end = end_x,
        anchor_ID = anchor_ID_x
    )],
    loops[, .(
        chr = chr_y,
        start = start_y,
        end = end_y,
        anchor_ID = anchor_ID_y
    )]
))

anchors_gr <- toGRanges(anchors)

acetyl_gr <- toGRanges(acetyl[, .SD, .SDcols = -c("width", "strand")])
acetyl_sig_gr <- toGRanges(acetyl_sig[, .SD, .SDcols = -c("width", "strand")])

hits <- as.data.table(findOverlaps(
    query = anchors_gr,
    subject = acetyl_gr
))

# sort hits by the loop type
hits <- merge(
    x = anchors[, .(.SD, .I)],        # include the index in the table for merging
    y = hits,
    by.x = "I",
    by.y = "queryHits"
)

# fix column names
# remove ".SD." prefix for column names
colnames(hits) <- gsub(".SD.", "anchor_", colnames(hits))
# remove "queryHits" column
hits[, I := NULL]

# merge in differential acetylation results
hits <- merge(
    x = hits,
    y = acetyl[, .(.SD, .I)],
    by.x = "subjectHits",
    by.y = "I",
    all.x = TRUE,
    all.y = FALSE
)

# fix column names
colnames(hits) <- gsub(".SD.", "peak_", colnames(hits))
# remove "subjectHits" column
hits[, subjectHits := NULL]

# repeat all of the above with only the differentially acetylated regions
hits_sig <- as.data.table(findOverlaps(
    query = anchors_gr,
    subject = acetyl_sig_gr
))

# sort hits by the loop type
hits_sig <- merge(
    x = anchors[, .(.SD, .I)],        # include the index in the table for merging
    y = hits_sig,
    by.x = "I",
    by.y = "queryHits"
)

# fix column names
# remove ".SD." prefix for column names
colnames(hits_sig) <- gsub(".SD.", "anchor_", colnames(hits_sig))
# remove "queryHits" column
hits_sig[, I := NULL]

# merge in differential acetylation results
hits_sig <- merge(
    x = hits_sig,
    y = acetyl_sig[, .(.SD, .I)],
    by.x = "subjectHits",
    by.y = "I",
    all.x = TRUE,
    all.y = FALSE
)

# fix column names
colnames(hits_sig) <- gsub(".SD.", "peak_", colnames(hits_sig))
# remove "subjectHits" column
hits_sig[, subjectHits := NULL]

loginfo("Performing calculations")


# ==============================================================================
# Save data
# ==============================================================================
loginfo("Saving data")
fwrite(
    hits[,
        .SD,
        .SDcols = c(
            "anchor_chr", "anchor_start", "anchor_end", "anchor_anchor_ID",
            "peak_chr", "peak_start", "peak_end", "peak_Conc_T2E", "peak_Conc_NonT2E", "peak_Fold",
            "peak_p.value", "peak_FDR"
        )
    ],
    "loops.intersected-acetyl.tsv",
    sep = "\t",
    col.names = TRUE
)

fwrite(
    hits_sig[,
        .SD,
        .SDcols = c(
            "anchor_chr", "anchor_start", "anchor_end", "anchor_anchor_ID",
            "peak_chr", "peak_start", "peak_end", "peak_Conc_T2E", "peak_Conc_NonT2E", "peak_Fold",
            "peak_p.value", "peak_FDR"
        )
    ],
    "loops.intersected-sig-acetyl.tsv",
    sep = "\t",
    col.names = TRUE
)
