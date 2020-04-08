# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("DESeq2"))
suppressMessages(library("GenomicRanges"))

GRAPH_DIR <- "Graphs"

# ==============================================================================
# Functions
# ==============================================================================
dt2gr <- function(dt) {
    '%ni%' <- Negate('%in%')
    non_pos_cols <- which(colnames(dt) %ni% c("chr", "start", "end"))
    gr <- GRanges(
        seqnames = dt$chr,
        ranges = IRanges(
            start = dt$start + 1,
            end = dt$end
        )
    )
    if (length(non_pos_cols) > 0) {
        mcols(gr) <- dt[, .SD, .SDcols = non_pos_cols]
    }
    return(gr)
}

split_comma_col <- function(v, f=identity) {
    # split into list
    splitv <- lapply(v, function(x) {strsplit(x, ",")[[1]]})
    # remove various non-informative characters (spaces, braces)
    splitv <- lapply(splitv, function(x) {gsub("[][ ]", "", x)})
    return(lapply(splitv, f))
}

# ==============================================================================
# Data
# ==============================================================================
# load sample metadata
metadata <- fread(
    "../../Data/External/LowC_Samples_Data_Available.tsv",
    sep = "\t",
    header = TRUE
)
metadata[, SampleID :=  paste0("PCa", get("Sample ID"))]
SAMPLES <- metadata$SampleID
metadata[, ChIP_file := paste0("../../Data/Processed/2019-05-03_PCa-H3K27ac-peaks/BAMs/Pca", get("Sample ID"), "_H3K27ac.sorted.dedup.bam")]
metadata[, Ctrl_file := paste0("../../Data/Processed/2019-05-03_PCa-H3K27ac-peaks/BAMs/Pca", get("Sample ID"), "_input.sorted.dedup.bam")]

# load TADs encompassing each breakpoint
tads <- fread(file.path(GRAPH_DIR, "sv-disruption-tests.TADs.tsv"), sep = "\t", header = TRUE)
tads_gr <- dt2gr(tads)

# load tests metadata
tests <- fread(file.path(GRAPH_DIR, "sv-disruption-tests.tsv"), sep = "\t", header = TRUE)
# convert comma-separated values into lists
tests$breakpoint_IDs <- split_comma_col(tests$breakpoint_IDs, as.numeric)
tests$mut_samples <- split_comma_col(tests$mut_samples)
tests$nonmut_samples <- split_comma_col(tests$nonmut_samples)
# assign placeholders for test data to be calculated
tests$z <- as.numeric(NA)
tests$p <- as.numeric(NA)
tests$padj <- as.numeric(NA)

# load induced regions (the rows of the count matrix)
regions <- fread(
    "TAD-induced-regions.bed",
    sep = "\t",
    header = FALSE,
    col.names = c("chr", "start", "end")
)
regions_gr <- dt2gr(regions)

# ==============================================================================
# Analysis
# ==============================================================================
# 1. Load ChIP-seq counts and calculate linear scale factors
# ----------------------------------------------------------
# load raw counts of H3K27ac pull down and input, then take differences and get total ChIP library sizes
count_matrix <- NULL
library_sizes <- sapply(SAMPLES, function(s) 0)
for (s in SAMPLES) {
    dt_chip <- fread(
        file.path("Acetylation", paste0(s, "_H3K27ac.induced-region-counts.bed")),
        sep = "\t",
        header = FALSE,
        col.names = c("chr", "start", "end", "count", "supported", "width", "frac_supported")
    )
    dt_input <- fread(
        file.path("Acetylation", paste0(s, "_input.induced-region-counts.bed")),
        sep = "\t",
        header = FALSE,
        col.names = c("chr", "start", "end", "count", "supported", "width", "frac_supported")
    )
    # calculate the linear scaling based on library size
    sample_library_size <- c(
        "chip" = dt_chip[, sum(count)],
        "input" = dt_input[, sum(count)]
    )
    scales <- c("chip" = 1, "input" = 1)
    if (sample_library_size["chip"] > sample_library_size["input"]) {
        scales["chip"] <- sample_library_size["input"] / sample_library_size["chip"]
        library_sizes[s] <- sample_library_size["input"]
    } else {
        scales["input"] <- sample_library_size["chip"] / sample_library_size["input"]
        library_sizes[s] <- sample_library_size["chip"]
    }

    # regions are in the same sorted order, so just taking difference
    dt_chip[, norm_count := pmax(
        1,
        ceiling(scales["chip"] * count - scales["input"] * dt_input$count)
    )]
    count_matrix <- cbind(count_matrix, dt_chip$norm_count)
}
colnames(count_matrix) <- SAMPLES

# 2. Calculate sample correlation based on scaled ChIP-seq counts
# ---------------------------------------------------------------
# calculcate correlations between acetylation in each sample
sample_corrs <- diag(ncol = length(SAMPLES), nrow = length(SAMPLES))
rownames(sample_corrs) <- SAMPLES
colnames(sample_corrs) <- SAMPLES
for (i in 1:(length(SAMPLES) - 1)) {
    for (j in (i + 1):length(SAMPLES)) {
        pairwise_cor <- cor(
            x = count_matrix[, i],
            y = count_matrix[, j],
            method = "spearman"
        )
        sample_corrs[i, j] <- pairwise_cor
        sample_corrs[j, i] <- pairwise_cor
    }
}

# 3. Perform differential acetylation testing for each group of TADs
# ------------------------------------------------------------------
# create design matrices for each set of mut-vs-nonmut (this depends on the test group)
# find all the unique combinations of foreground and background samples (mut vs nonmut)
# this will make creating the design matrices easier
test_combos <- lapply(
    1:tests[, .N],
    function(i) {
        return(list(
            "mut" = tests[i, unlist(mut_samples)],
            "nonmut" = tests[i, unlist(nonmut_samples)]
        ))
    }
)
all_comparisons <- unique(test_combos)

# map all the unique comparisons back to the corresponding test IDs
for (i in 1:length(all_comparisons)) {
    comp_mut <- all_comparisons[[i]]$mut
    comp_nonmut <- all_comparisons[[i]]$nonmut
    for (j in 1:tests[, .N]) {
        test_mut <- tests[j, unlist(mut_samples)]
        test_nonmut <- tests[j, unlist(nonmut_samples)]
        # only add this test_ID if the foreground and background mutations match
        # and this region is being tested in the first place
        if (identical(test_mut, comp_mut) && identical(test_nonmut, comp_nonmut)) {
            all_comparisons[[i]]$test_IDs <- c(all_comparisons[[i]]$test_IDs, tests[j, test_ID])
        }
    }
}

# perform differential analysis for each test
for (i in 1:length(all_comparisons)) {
    cat(i, "of", length(all_comparisons), "\n")
    mut_samples <- all_comparisons[[i]]$mut
    nonmut_samples <- all_comparisons[[i]]$nonmut
    testing_samples <- c(mut_samples, nonmut_samples)
    test_IDs <- all_comparisons[[i]]$test_IDs
    # create metadata table for samples
    meta <- data.table(
        Sample = testing_samples,
        Mutated = rep(
            c("Yes", "No"),
            c(length(mut_samples), length(nonmut_samples))
        ),
        Size = library_sizes[testing_samples]
    )
    meta[, Mutated := factor(Mutated, levels = c("No", "Yes"))]
    # create DESeqDataSet object
    dds <- DESeqDataSetFromMatrix(
        countData = count_matrix[, testing_samples],
        colData = meta,
        design = ~ Mutated,
        rowRanges = regions_gr
    )
    # filter out regions with few reads in the region
    cat(">>filtering regions\n")
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep, ]
    cat(">>performing normalization\n")
    # perform steps indivudally, instead of having them wrapped in `DESeq` function
    # as per DiffBind documentation and other ChIP-seq analysis tools, linearly scale according to the sample with
    # the smallest library
    sizeFactors(dds) <- library_sizes[testing_samples] / min(library_sizes[testing_samples])
    dds <- estimateDispersions(dds, fitType="local")
    dds <- nbinomWaldTest(dds)
    cat(">>extracting results\n")
    res <- results(dds)
    res_dt <- data.table(
        "chr" = as.character(seqnames(granges(dds))),
        "start" = start(granges(dds)) - 1,
        "end" = end(granges(dds)),
        "baseMean" = res$baseMean,
        "log2FoldChange" = res$log2FoldChange,
        "lfcSE" = res$lfcSE,
        "stat" = res$stat,
        "p" = res$pvalue,
        "padj" = res$padj
    )

    # the above is done on all regions to perform proper normalization
    # we're only interesting in testing the TADs containing the breakpoints
    # we need to find the indices of all the regions that overlap the TADs of interest
    # for this particular contrast

    # find all the tests with the contrast of interest
    # enforce positions as integers to avoid writing to the file with scientific notation
    res_dt[, start := as.integer(start)]
    res_dt[, end := as.integer(end)]
    # save test results
    for (ti in test_IDs) {
        fwrite(
            res_dt,
            paste0("Acetylation/Tests/test_", ti, ".results.tsv"),
            sep = "\t",
            col.names = TRUE
        )
    }

    # find all the TADs related to these tests
    tads_of_interest <- tads_gr[tads_gr$test_ID %in% test_IDs]

    # find all the induced regions that map to these TADs
    regions_of_interest <- as.data.table(findOverlaps(regions_gr, tads_of_interest))
    regions_of_interest[, test_ID := tads_of_interest[subjectHits]$test_ID]
    # get differential analysis results for these regions
    regions_of_interest[, p := res_dt[queryHits, p]]
    # calculate region weights for the TAD it is induced from
    regions_of_interest[, weight := sqrt(
        width(GenomicRanges::pintersect(
            regions_gr[queryHits],
            tads_of_interest[subjectHits]
        )) / width(tads_of_interest[subjectHits])
    )]
    # calculate Stouffer's Z with the weights calculated above
    stouffer_z <- regions_of_interest[,
        #!is.na(p) & !is.na(padj),
        sum(weight * qnorm(1 - p), na.rm = TRUE) / sqrt(sum(weight ^ 2)),
        by = test_ID
    ]
    colnames(stouffer_z)[2] <- "z"
    # for any test that has all NAs, make sure it's listed as such (because the sum will be 0)
    tests_with_all_nas <- regions_of_interest[, all(is.na(p)), by = "test_ID"][V1 == TRUE, test_ID]
    stouffer_z[test_ID %in% tests_with_all_nas, z := NA]
    # calculate the associated p-value
    stouffer_z[, p := 2 * pnorm(-abs(z))]
    # because we're only testing these regions, not the others, perform p-value adjustment here
    stouffer_z[, padj := p.adjust(p, method = "fdr", n = .N)]

    # save to table for later
    for (ti in test_IDs) {
        tests[test_ID == ti, z := stouffer_z[test_ID == ti]$z]
        tests[test_ID == ti, p := stouffer_z[test_ID == ti]$p]
        tests[test_ID == ti, padj := stouffer_z[test_ID == ti]$padj]
        tests[test_ID == ti, chr := as.character(seqnames(tads_of_interest[tads_of_interest$test_ID == ti]))]
        tests[test_ID == ti, start := start(tads_of_interest[tads_of_interest$test_ID == ti]) - 1]
        tests[test_ID == ti, end := end(tads_of_interest[tads_of_interest$test_ID == ti])]
    }
}

# ==============================================================================
# Save data
# ==============================================================================
# convert lists back to comma-separated columns
tests_backup <- copy(tests)
tests$breakpoint_IDs <- unlist(lapply(tests$breakpoint_IDs, paste, collapse = ","))
tests$mut_samples <- unlist(lapply(tests$mut_samples, paste, collapse = ","))
tests$nonmut_samples <- unlist(lapply(tests$nonmut_samples, paste, collapse = ","))

# enforce positions as integers to avoid writing to the file with scientific notation
tests[, start := as.integer(start)]
tests[, end := as.integer(end)]

# write the tables with columns in a particular order
fwrite(
    tests[, .SD, .SDcols = c("test_ID", "z", "p", "padj")],
    "Graphs/sv-disruption-tests.acetylation.tsv",
    sep = "\t",
    col.names = TRUE
)
fwrite(
    sample_corrs,
    "Statistics/acetylation-correlation.tsv",
    sep = "\t",
    col.names = TRUE
)
