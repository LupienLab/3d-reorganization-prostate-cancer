# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
suppressMessages(library("GenomicRanges"))
suppressMessages(library("regioneR"))

PEAK_DIR <- file.path("..", "..", "Data", "Processed", "2019-05-03_PCa-H3K27ac-peaks", "Peaks")
ACETYL_DIR <- file.path("Acetylation", "Tests")
EXPRS_DIR <- file.path("..", "2020-06-18_sv-disruption-expression", "sleuth")

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
gencode_gr <- toGRanges(gencode, genome = "hg38")
# fix improper strand parsing
strand(gencode_gr) <- gencode[, strand]
mcols(gencode_gr)$strand <- NULL


# ==============================================================================
# Analysis
# ==============================================================================
# calculate weighted-mean gene expression fold change for each gene in each test
delta_exprs <- exprs[, .(exprs = sum(mean_obs * b) / sum(mean_obs)), keyby = c("test_ID", "ens_gene", "gene_name")]

# map peaks to genes they overlap
hits <- as.data.table(findOverlaps(acetyl_gr, gencode_gr))
matched_acetyl_gr <- acetyl_gr[hits$queryHits]
strand(matched_acetyl_gr) <- strand(gencode_gr[hits$subjectHits])
matched_acetyl_gr$gene_id <- gencode_gr[hits$subjectHits]$gene_id
matched_acetyl_gr$gene_name <- gencode_gr[hits$subjectHits]$gene_name

matched_acetyl <- as.data.table(matched_acetyl_gr)

# calculate weighted-mean acetylation fold change for each gene in each test
delta_acetyl <- matched_acetyl[, .(acetyl = sum(Conc_Nonmutated * Fold) / sum(Conc_Nonmutated)), keyby = c("test_ID", "gene_id", "gene_name")]

# merge differential results into the same table
delta <- merge(
    x = delta_exprs,
    y = delta_acetyl,
    by.x = c("test_ID", "ens_gene", "gene_name"),
    by.y = c("test_ID", "gene_id", "gene_name"),
    all = FALSE
)

# ==============================================================================
# Save data
# ==============================================================================
fwrite(delta, "exprs-acetyl.tsv", sep = "\t", col.names = TRUE)
