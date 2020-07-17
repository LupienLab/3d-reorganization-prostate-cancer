# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
suppressMessages(library("GenomicRanges"))
suppressMessages(library("regioneR"))
suppressMessages(library("DiffBind"))

PEAK_DIR <- file.path("..", "..", "Data", "Processed", "2019-05-03_PCa-H3K27ac-peaks", "Peaks")
ACETYL_DIR <- file.path("Acetylation", "Tests")
EXPRS_DIR <- file.path("..", "2020-06-18_sv-disruption-expression", "sleuth")

QVAL_THRESH <- 0.05

# ==============================================================================
# Data
# ==============================================================================
# load sample metadata
metadata <- fread(
    "config.tsv",
    sep = "\t",
    header = TRUE
)
low_qual_samples <- metadata[Include == "No"]
metadata <- metadata[Include == "Yes"]
SAMPLES <- metadata$SampleID

# test IDs for differential acetylation data
test_IDs <- sort(sapply(
    list.files(ACETYL_DIR, "local.tsv"),
    function(fn) {
        as.numeric(gsub("\\.local\\.tsv", "", gsub("test_", "", fn)))
    }
))

# differential acetylation results
acetyl <- rbindlist(lapply(
    test_IDs,
    function(tid) {
        dt <- fread(
            file.path(ACETYL_DIR, paste0("test_", tid, ".local.tsv")),
            sep = "\t",
            header = TRUE,
            drop = c("strand", "width")
        )
        dt[, test_ID := tid]
        return(dt)
    }
))
acetyl_gr <- toGRanges(acetyl, genome = "hg38")
strand(acetyl_gr) <- "+"
acetyl_gr_copy <- copy(acetyl_gr)
strand(acetyl_gr_copy) <- "-"
acetyl_gr <- c(acetyl_gr, acetyl_gr_copy)
rm(acetyl_gr_copy)

# differential expression results
exprs <- rbindlist(lapply(
    test_IDs,
    function(tid) {
        dt <- fread(
            file.path(EXPRS_DIR, paste0("test_", tid, ".transcripts.tested.tsv")),
            sep = "\t",
            header = TRUE
        )
        dt[, test_ID := tid]
        # multiple test correction since we're only looking at these transcripts and genes in the first place
        dt[, qval := p.adjust(pval, method = "fdr")]
        return(dt)
    }
))
exprs_genes <- rbindlist(lapply(
    test_IDs,
    function(tid) {
        dt <- fread(
            file.path(EXPRS_DIR, paste0("test_", tid, ".genes.tested.tsv")),
            sep = "\t",
            header = TRUE
        )
        dt[, test_ID := tid]
        # multiple test correction since we're only looking at these transcripts and genes in the first place
        dt[, qval := p.adjust(pval, method = "fdr")]
        return(dt)
    }
))

# load GENCODE reference
gencode <- fread(
    file.path("..", "..", "Data", "External", "GENCODE", "gencode.v33.all-genes.bed"),
    sep = "\t",
    header = FALSE,
    col.names = c("chr", "start", "end", "strand", "gene_id", "gene_name")
)

# shift gene locations to include promoter regions
gencode[strand == "+", start := start - 5000]
gencode[strand == "-", end := end + 5000]

# convert to GRanges
gencode_gr <- toGRanges(gencode, genome = "hg38")
# fix improper strand parsing
strand(gencode_gr) <- gencode[, strand]
mcols(gencode_gr)$strand <- NULL

# load DiffBind data
dba_all <- readRDS("dba_all.rds")

# ==============================================================================
# Analysis
# ==============================================================================
# calculate weighted-mean gene expression fold change for each gene in each test
# `b` column from sleuth is natural log, convert to log2 before saving
delta_exprs <- exprs[,
    .(log2exprs_fc = sum(mean_obs * b * log2(exp(1))) / sum(mean_obs)),
    keyby = "ens_gene"
]
delta_exprs <- merge(
    x = delta_exprs,
    y = exprs_genes[, .SD, .SDcols = c("target_id", "gene_name", "pval", "qval", "test_ID")],
    by.x = "ens_gene",
    by.y = "target_id"
)
colnames(delta_exprs)[colnames(delta_exprs) == "pval"] <- "exprs_pval"
colnames(delta_exprs)[colnames(delta_exprs) == "qval"] <- "exprs_qval"

# map peaks to genes they overlap
near_hits <- nearest(acetyl_gr, gencode_gr)
acetyl_gr$nearest_gene_id <- gencode_gr[near_hits]$gene_id

hits <- as.data.table(findOverlaps(acetyl_gr, gencode_gr))
matched_acetyl_gr <- acetyl_gr[hits$queryHits]
matched_acetyl_gr$gene_id <- gencode_gr[hits$subjectHits]$gene_id
matched_acetyl_gr$gene_name <- gencode_gr[hits$subjectHits]$gene_name
matched_acetyl <- as.data.table(matched_acetyl_gr)

# calculate weighted-mean acetylation fold change for each gene in each test
delta_acetyl <- matched_acetyl[,
    .(
        log2acetyl_fc = sum(Conc_Nonmutated * Fold) / sum(Conc_Nonmutated),
        acetyl_pval = pchisq((-2 * sum(log(p.value))), df = 2 * .N, lower.tail = FALSE),
        # Fisher's method for combination, followed by FDR correction
        acetyl_qval = p.adjust(pchisq((-2 * sum(log(p.value))), df = 2 * .N, lower.tail = FALSE), method = "fdr")
    ),
    keyby = c("gene_id", "gene_name")
    # keyby = c("gene_id", "gene_name", "p.value", "FDR")
]


# merge differential results into the same table
delta <- merge(
    x = delta_exprs,
    y = delta_acetyl,
    by.x = c("ens_gene", "gene_name"),
    by.y = c("gene_id", "gene_name"),
    all = FALSE
)
delta[, sig := ((exprs_qval <  QVAL_THRESH) & (acetyl_qval <  QVAL_THRESH))]

# ==============================================================================
# Save data
# ==============================================================================
fwrite(delta, "exprs-acetyl.tsv", sep = "\t", col.names = TRUE)
