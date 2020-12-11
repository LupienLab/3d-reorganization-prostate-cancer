# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("logging"))

loginfo("Loading packages")
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
suppressMessages(library("regioneR"))
suppressMessages(library("BSgenome.Hsapiens.UCSC.hg38.masked"))
source("../2020-02-19_chromoplexy/plotting-helper.R")

N_PERMS <- 1000


# ==============================================================================
# Data
# ==============================================================================
loginfo("Loading info")

# load CTCF binding sites in prostate cancer cell lines
CELL_LINES <- c("22Rv1", "LNCaP", "C4-2B", "LNCaP", "VCaP")
pca_ctcf <- rbindlist(lapply(
    CELL_LINES,
    function(cl) {
        dt <- fread(
            paste0("../../Data/External/ENCODE_ChIP/", cl, "_CTCF_Optimal.canonical.sorted.bed"),
            sep = "\t",
            header = FALSE,
            select = c(1:3, 5),
            col.names = c("chr", "start", "end", "score")
        )
        dt[, Cell_Line := cl]
        return(dt)
    }
))
ctcf <- GRangesList(lapply(
    CELL_LINES,
    function(cl) {
        ctcf_bs <- pca_ctcf[Cell_Line == cl]
        return(GRanges(
            seqnames = ctcf_bs$chr,
            ranges = IRanges(
                start = ctcf_bs$start + 1,
                end = ctcf_bs$end
            )
        ))
    }
))
names(ctcf) <- CELL_LINES

# load DMRs
pca_dmrs <- fread("../../Data/External/Zhao_DNAme/TableS8_DMRs_PCa-vs-benign.tsv", sep = "\t", header = TRUE)
# only look at DMRs that are hypermethylated in prostate cancer, compared to benigns
pca_dmrs <- pca_dmrs[get("difference in methylation") > 0]
dmrs <- GRanges(
    seqnames = pca_dmrs$chr,
    ranges = IRanges(
        start = pca_dmrs$start + 1,
        end = pca_dmrs$end
    ),
    delta_beta = pca_dmrs$`difference in methylation`,
    area = pca_dmrs$areaStat
)


# ==============================================================================
# Analysis
# ==============================================================================
loginfo("Loading genome annotation")

# get mask of regions to exclude from permutations
hg38 <- getGenomeAndMask("hg38")
hg38_canonical <- filterChromosomes(hg38$genome, organism="hg", chr.type = "canonical")
hg38_mask <- filterChromosomes(hg38$mask, organism="hg", chr.type = "canonical")

loginfo("Performing permutation tests")
# permutation test for number of overlaps
perm_test <- lapply(
    CELL_LINES,
    function(cl) {
        overlapPermTest(
            A = dmrs,
            B = ctcf[[cl]],
            genome = hg38_canonical,
            ntimes = N_PERMS,
            alternative = "greater",
            mask = hg38_mask
        )
    }
)
names(perm_test) <- CELL_LINES

perm_test_data <- data.table(Cell_Line = CELL_LINES)
perm_test_data[, observed := apply(.SD, 1, function(r) perm_test[[r["Cell_Line"]]]$numOverlaps$observed)]
perm_test_data[, p := apply(.SD, 1, function(r) perm_test[[r["Cell_Line"]]]$numOverlaps$pval)]
perm_test_data[, FDR := p.adjust(p, method = "fdr")]

loginfo("Local adjustment for sensitivity")
# local positioning adjustment to assess sensitivity to exact placement
lz <- lapply(
    CELL_LINES,
    function(cl) {
        localZScore(
            pt = perm_test[[cl]],
            A = dmrs,
            B = ctcf[[cl]],
            ntimes = N_PERMS,
            window = 5000,
            step = 10
        )
    }
)
names(lz) <- CELL_LINES

# local enrichment sensitivity data
lz_dt <- rbindlist(lapply(
    CELL_LINES,
    function(cl) {
        data.table(
            Cell_Line = cl,
            Shift = lz[[cl]]$numOverlaps$shifts,
            z = lz[[cl]]$numOverlaps$shifted.z.scores
        )
    }
))

# ==============================================================================
# Save data
# ==============================================================================
loginfo("Saving permutation data")

fwrite(
    perm_test_data,
    "Overlaps/pca-hyper-dmrs.permutations.tsv",
    sep = "\t",
    col.names = TRUE
)

fwrite(
    lz_dt,
    "Overlaps/pca-hyper-dmrs.local-dependency.tsv",
    sep = "\t",
    col.names = TRUE
)

saveRDS(perm_test, "Overlaps/pca-hyper-dmrs.permutation-tests.rds")
saveRDS(lz, "Overlaps/pca-hyper-dmrs.local-dependency.rds")
