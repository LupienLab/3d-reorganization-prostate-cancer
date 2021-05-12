# Perform hypothesis testing to compare the SV breakpoints detected by WGS and Hi-C
# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
suppressMessages(library("GenomicRanges"))

RES_DIR = path.join("..", "..", "results", "2020-02-20_wgs-hic-sv-comparison")

# ==============================================================================
# Functions
# ==============================================================================
dt2gr <- function(dt) {
    possible_chr_col_names <- c("chr", "seqnames")
    possible_start_col_names <- c("start", "from", "x1")
    possible_end_col_names <- c("end", "to", "x2")
    chr_col <- possible_chr_col_names[possible_chr_col_names %in% colnames(dt)][1]
    start_col <- possible_start_col_names[possible_start_col_names %in% colnames(dt)][1]
    end_col <- possible_end_col_names[possible_end_col_names %in% colnames(dt)][1]
    gr <- GRanges(
        seqnames = dt[, get(chr_col)],
        ranges = IRanges(
            start = dt[, get(start_col) + 1],
            end = dt[, get(end_col)],
        )
    )
    mcols(gr) <- dt[,
        .SD,
        .SDcols = setdiff(colnames(dt), c(chr_col, start_col, end_col))
    ]
    return(gr)
}

# ==============================================================================
# Data
# ==============================================================================
# load metadata
metadata <- fread(
    file.path("..", "..", "data", "External", "LowC_Samples_Data_Available.tsv"),
    sep = "\t",
    header = TRUE
)
metadata <- metadata[Include == "Yes"]
metadata[, SampleID := paste0("PCa", get("Sample ID"))]
SAMPLES <- metadata[, SampleID]

# load Hi-C-detected breakpoints
hic <- fread(
    file.path("..", "..", "results", "2020-02-19_chromoplexy", "Graphs", "sv-breakpoints.tsv"),
    header = TRUE,
    sep = "\t"
)
hic_gr <- dt2gr(hic)

# load Hi-C-detected breakpoint pairs
hic_pairs <- fread(
    file.path("..", "..", "results", "2020-02-19_chromoplexy", "Graphs", "sv-breakpoints.paired.tsv"),
    header = TRUE,
    sep = "\t"
)

# load WGS-detected breakpoints
wgs <- rbindlist(lapply(
    SAMPLES,
    function(s) {
        dt <- fread(
            file.path(
                "..", "..", "data", "External", "CPC-GENE", "structural-variant-vcfs",
                paste0(s, ".hg38.sorted.bed")
            ),
            header = FALSE,
            sep = "\t",
            col.names = c("chr", "start", "end", "SV_ID")
        )
        dt[, SampleID := s]
        return(dt)
    }
))
wgs[, breakpoint_ID := paste(SampleID, SV_ID, .I, sep = "_")]
wgs_gr <- dt2gr(wgs)

# load ATAC-seq TCGA peaks for all tissue types
atac <- rbindlist(lapply(
    list.files(
        file.path("..", "..", "data", "External", "TCGA"),
        pattern = "peakCalls.txt",
        full.names = TRUE
    ),
    function(f) {
        dt <- fread(f, header = TRUE, sep = "\t")
        dt[, Tissue := substr(
            basename(f),
            1,
            regexpr("_", basename(f)
        )[[1]] - 1)]
    }
))
ATAC_TISSUES <- atac[, unique(Tissue)]
atac_gr <- GRangesList(lapply(
    ATAC_TISSUES,
    function(t) dt2gr(atac[Tissue == t])
))
names(atac_gr) <- ATAC_TISSUES

# ==============================================================================
# Analysis
# ==============================================================================

# 1. Compare Hi-C and WGS breakpoints by their chromatin accessibility
# --------------------------------------
# get median size of breakpoints
EXT_SIZE <- median(width(hic_gr))

wgs_gr_extended <- copy(wgs_gr)
start(wgs_gr_extended) <- start(wgs_gr) - floor(EXT_SIZE / 2)
end(wgs_gr_extended) <- end(wgs_gr) + floor(EXT_SIZE / 2)

# find overlaps between breakpoints and ATAC peaks

# try it for all tissues to see if there is any difference due to the
# ATAC peak locations in different tissues
hic_accessible_idx_all_tissues <- sapply(
    ATAC_TISSUES,
    function(t) {
        length(unique(queryHits(findOverlaps(hic_gr, atac_gr[[t]]))))
    },
    USE.NAMES = TRUE
)
wgs_accessible_idx_all_tissues <- sapply(
    ATAC_TISSUES,
    function(t) {
        length(unique(queryHits(findOverlaps(wgs_gr, atac_gr[[t]]))))
    },
    USE.NAMES = TRUE
)

hic_accessible_idx <- unique(queryHits(findOverlaps(hic_gr, atac_gr$PRAD)))
wgs_accessible_idx <- unique(queryHits(findOverlaps(wgs_gr, atac_gr$PRAD)))
wgs_extended_accessible_idx <- unique(
    queryHits(findOverlaps(wgs_gr_extended, atac_gr$PRAD))
)

hic[, Accessible := FALSE]
hic[hic_accessible_idx, Accessible := TRUE]
wgs[, Accessible := FALSE]
wgs_extended <- copy(wgs)
wgs[wgs_accessible_idx, Accessible := TRUE]
wgs_extended[wgs_extended_accessible_idx, Accessible := TRUE]

# convert into multi-dimensional array that works with CMH test
accessibility_counts <- unlist(lapply(
    SAMPLES,
    function(s) c(
        wgs[SampleID == s & Accessible == TRUE, .N],
        wgs[SampleID == s & Accessible == FALSE, .N],
        hic[SampleID == s & Accessible == TRUE, .N],
        hic[SampleID == s & Accessible == FALSE, .N]
    )
))
accessibility_array <- array(
    accessibility_counts,
    dim = c(2, 2, 12),
    dimnames = list(
        Accessible = c(TRUE, FALSE),
        Method = c("WGS", "Hi-C"),
        SampleID = SAMPLES
    )
)

# do the same thing with extended WGS breakpoints
accessibility_extended_counts <- unlist(lapply(
    SAMPLES,
    function(s) c(
        wgs_extended[SampleID == s & Accessible == TRUE, .N],
        wgs_extended[SampleID == s & Accessible == FALSE, .N],
        hic[SampleID == s & Accessible == TRUE, .N],
        hic[SampleID == s & Accessible == FALSE, .N]
    )
))
accessibility_extended_array <- array(
    accessibility_extended_counts,
    dim = c(2, 2, 12),
    dimnames = list(
        Accessible = c(TRUE, FALSE),
        Method = c("WGS", "Hi-C"),
        SampleID = SAMPLES
    )
)

# 2. Compare Hi-C and WGS breakpoint pair types
# --------------------------------------
ALL_SV_TYPES <- c("INV", "DUP", "DEL", "BND", "UNKNOWN")

wgs_pairs <- wgs[, .(SV_Type = substr(unique(SV_ID), 1, 3)), by = "SampleID"]

sv_types <- as.data.table(expand.grid(
    Method = c("WGS", "Hi-C"),
    SV_Type = ALL_SV_TYPES,
    SampleID = SAMPLES
))
sv_types[, N := apply(.SD, 1, function(r) {
    ifelse(
        r["Method"] == "WGS",
        wgs_pairs[SampleID == r["SampleID"] & SV_Type == r["SV_Type"], .N],
        hic_pairs[SampleID == r["SampleID"] & sv_type == r["SV_Type"], .N]
    ) 
})]

bnd_counts <- unlist(lapply(
    SAMPLES,
    function(s) unlist(lapply(
        c("WGS", "Hi-C"),
        function(m) c(
            sv_types[SampleID == s & Method == m & SV_Type == "BND", sum(N)],
            sv_types[SampleID == s & Method == m & SV_Type != "BND", sum(N)]
        )
    ))
))
bnd_array <- array(
    bnd_counts,
    dim = c(2, 2, 12),
    dimnames = list(
        BND = c(TRUE, FALSE),
        Method = c("WGS", "Hi-C"),
        SampleID = SAMPLES
    )
)


# 3. Perform hypothesis tests for each of the above covariates
# --------------------------------------
htests <- list(
    Accessibility = mantelhaen.test(accessibility_array, alternative = "less"),
    Accessibility_Extended = mantelhaen.test(
        accessibility_extended_array,
        alternative = "less"
    ),
    SV_Type = mantelhaen.test(bnd_array, alternative = "less")
)

# ==============================================================================
# Save data
# ==============================================================================
htests_dt <- rbindlist(lapply(
    names(htests),
    function(n) data.table(
        Test = n,
        Variable = names(htests[[n]]),
        Value = htests[[n]]
    )
))
fwrite(
    htests_dt,
    path.join(RES_DIR, "sv-annotation-comparisons.tsv"),
    sep = "\t",
    col.names = TRUE
)

fwrite(
    sv_types,
    path.join(RES_DIR, "sv-types-counted.tsv"),
    sep = "\t",
    col.names = TRUE
)
