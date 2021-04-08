# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
suppressMessages(library("regioneR"))
suppressMessages(library("BSgenome.Hsapiens.UCSC.hg38.masked"))


# ==============================================================================
# Data
# ==============================================================================
# load SV breakpoints
breakpoints <- fread(
    file.path("enrichment", "sv-breakpoints.bed"),
    sep = "\t",
    header = FALSE,
    col.names = c("chr", "start", "end")
)
# add +1 for start offset
breakpoints[, start := start + 1]
bp_gr <- toGRanges(breakpoints, genome = "hg38")

# load H3K27ac peaks
h3k27ac <- fread(
    file.path("..", "..", "Data", "Processed", "2019-05-03_PCa-H3K27ac-peaks", "Peaks", "catalogue-peaks.bed"),
    sep = "\t",
    header = FALSE,
    col.names = c("chr", "start", "end")
)
h3k27ac_gr <- toGRanges(h3k27ac, genome = "hg38")

# load ENCODE blacklist
blacklist <- fread(
    file.path("..", "..", "Data", "External", "ENCODE_ChIP", "ENCODE-blacklist.bed"),
    sep = "\t",
    header = FALSE,
    col.names = c("chr", "start", "end")
)
blacklist_gr <- toGRanges(blacklist, genome = "hg38")

# ==============================================================================
# Analysis
# ==============================================================================
# get mask of regions to exclude from permutations
hg38 <- getGenomeAndMask("hg38")
hg38_canonical <- filterChromosomes(hg38$genome, organism="hg", chr.type = "canonical")
hg38_mask <- filterChromosomes(hg38$mask, organism="hg", chr.type = "canonical")
hg38_mask <- mergeRegions(hg38_mask, blacklist_gr)

# permutation test for number of overlaps between SV breakpoints and H3K27ac peaks
perm_test <- overlapPermTest(
    A = bp_gr,
    B = h3k27ac_gr,
    genome = hg38_canonical,
    ntimes = 100,
    alternative = "greater",
    mask = hg38_mask
)

# local positioning adjustment to assess sensitivity to exact placement
lz <- localZScore(
    pt = perm_test,
    A = bp_gr,
    B = h3k27ac_gr,
    ntimes = 100,
    window = 10000,
    step = 1000
)

perm_test_data <- data.table(
    N_Iterations = perm_test$numOverlaps$ntimes,
    observed_intersections = perm_test$numOverlaps$observed,
    mean_permuted_intersections = mean(perm_test$numOverlaps$permuted),
    sd_permuted_intersections = sd(perm_test$numOverlaps$permuted),
    z = perm_test$numOverlaps$zscore,
    p = perm_test$numOverlaps$pval
)

# ==============================================================================
# Plots
# ==============================================================================
png(
    file.path("Plots", "permutation.png"),
    width = 20,
    height = 12,
    units = "cm",
    res = 300
)
plot(perm_test)
dev.off()
pdf(
    file.path("Plots", "permutation.pdf"),
    width = 20,
    height = 12
)
plot(perm_test)
dev.off()

png(
    file.path("Plots", "local-z.png"),
    width = 20,
    height = 12,
    units = "cm",
    res = 300
)
plot(lz)
dev.off()
pdf(
    file.path("Plots", "local-z.pdf"),
    width = 20,
    height = 12
)
plot(lz)
dev.off()

# ==============================================================================
# Save data
# ==============================================================================

fwrite(
    perm_test_data,
    file.path("enrichment", "sv-breakpoints.H3K27ac-catalogue.permutations.tsv"),
    sep = "\t",
    col.names = TRUE
)

saveRDS(perm_test, file.path("enrichment", "permutation-tests.rds"))
saveRDS(lz, file.path("enrichment", "local-dependency.rds"))
