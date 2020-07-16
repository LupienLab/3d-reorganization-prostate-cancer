# ==============================================================================
# Environment
# ==============================================================================
cat("Loading packages\n")
suppressMessages(library("data.table"))
suppressMessages(library("DESeq2"))
suppressMessages(library("DiffBind"))
suppressMessages(library("GenomicRanges"))

BAM_DIR <- file.path("..", "..", "Data", "Processed", "2019-05-03_PCa-H3K27ac-peaks", "BAMs")
PEAK_DIR <- file.path("..", "..", "Data", "Processed", "2019-05-03_PCa-H3K27ac-peaks", "Peaks")
TEST_DIR <- file.path("Acetylation", "T2E")


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
metadata[, `:=`(
    bamReads = file.path(BAM_DIR, paste0(gsub("PCa", "Pca", SampleID), "_H3K27ac.sorted.dedup.bam")),
    bamControl = file.path(BAM_DIR, paste0(gsub("PCa", "Pca", SampleID), "_input.sorted.dedup.bam")),
    Peaks = file.path(PEAK_DIR, paste0(SampleID, "_peaks.filtered.sorted.merged.narrowPeak")),
    PeakCaller = "narrow"
)]


# ==============================================================================
# Analysis
# ==============================================================================
# initialize DBA object
cat("Loading metadata\n")
dba_all <- dba(sampleSheet = metadata)
# quantify reads in each sample
cat("Quantifying reads\n")
dba_all <- dba.count(dba_all, summits=250)
saveRDS(dba_all, "dba_all.rds")

# establish contrast for differential analysis
cat("Differential analysis\n")
dba_comp <- dba.contrast(
    dba_all,
    group1 = metadata[, T2E_Status == "Yes"],
    group2 = metadata[, T2E_Status == "No"],
    name1 = "T2E",
    name2 = "NonT2E"
)

# perform differential analysis
dba_comp <- dba.analyze(dba_comp, bReduceObjects = FALSE)
saveRDS(dba_comp, file.path(TEST_DIR, "dba_comp.rds"))

# extract differntially acetylated regions
comp_regions <- dba.report(dba_comp, th = 1)

# convert results to data.table
comp_dt <- as.data.table(comp_regions)
comp_dt[, start := start - 1]
colnames(comp_dt)[1] <- "chr"


# ==============================================================================
# Save data
# ==============================================================================
cat("Saving data\n")
# write the full set of results
fwrite(
    comp_dt,
    file.path(TEST_DIR, "t2e.all.tsv"),
    sep = "\t",
    col.names = TRUE
)

# get loci only in the local region (i.e. the test TAD)
local_region <- list(
    chr = "chr21",
    start = 30000000,
    end = 50000000
)
local_comp <- comp_dt[(chr == local_region$chr) & (start <= local_region$end) & (end >= local_region$start)]

fwrite(
    local_comp,
    file.path(TEST_DIR, "t2e.local.tsv"),
    sep = "\t",
    col.names = TRUE
)
