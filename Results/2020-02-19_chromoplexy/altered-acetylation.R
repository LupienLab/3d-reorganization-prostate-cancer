# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("DESeq2"))
suppressMessages(library("GenomicRanges"))

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
SAMPLES <- paste0("PCa", metadata[, get("Sample ID")])
metadata[, ChIP_file := paste0("../../Data/Processed/2019-05-03_PCa-H3K27ac-peaks/BAMs/Pca", get("Sample ID"), "_H3K27ac.sorted.dedup.bam")]
metadata[, Ctrl_file := paste0("../../Data/Processed/2019-05-03_PCa-H3K27ac-peaks/BAMs/Pca", get("Sample ID"), "_input.sorted.dedup.bam")]

# load TADs encompassing each breakpoint
tads <- fread("sv-disruption-tests.TADs.tsv", sep = "\t", header = TRUE)
tads$test_index <- split_comma_col(tads$test_index, as.numeric)
tads[, any_altered_TAD := grepl("True", altered_TAD)]
tads$altered_TAD <- split_comma_col(tads$altered_TAD, function(x) grepl("True", x))
tads_gr <- dt2gr(tads)

# load tests metadata
tests <- fread("sv-disruption-tests.tsv", sep = "\t", header = TRUE)
# convert comma-separated values into lists
tests$breakpoint_indices <- split_comma_col(tests$breakpoint_indices, as.numeric)
tests$mutated_in <- split_comma_col(tests$mutated_in)
# assign placeholders for test data to be calculated
tests$z <- as.numeric(NA)
tests$p <- as.numeric(NA)
tests$padj <- as.numeric(NA)
tests$chr <- as.character(NA)
tests$start <- as.numeric(NA)
tests$end <- as.numeric(NA)

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
# create design matrices for each set of mut-vs-nonmut (this depends on the breakpoint being considered)
all_comparisons <- unique(tests$mutated_in)

# perform differential analysis for each test
for (i in 1:length(all_comparisons)) {
    cat(i, "of", length(all_comparisons), "\n")
    mut_samples <- all_comparisons[[i]]
    # create metadata table for samples
    meta <- data.table(Sample = SAMPLES, Mutated = "No", Size = library_sizes)
    meta[Sample %in% mut_samples, Mutated := "Yes"]
    meta[, Mutated := factor(Mutated, levels = c("No", "Yes"))]
    # create DESeqDataSet object
    dds <- DESeqDataSetFromMatrix(
        countData = count_matrix,
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
    sizeFactors(dds) <- library_sizes / min(library_sizes)
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
    # (i.e. mut_samples is equal to the `mutated_in` column)
    test_idx_of_interest <- tests[
        which(apply(tests, 1, function(r) identical(r$mutated_in, mut_samples))),
        test_index
    ]
    # enforce positions as integers to avoid writing to the file with scientific notation
    res_dt[, start := as.integer(start)]
    res_dt[, end := as.integer(end)]
    # save test results
    for (ti in test_idx_of_interest) {
        fwrite(
            res_dt,
            paste0("Acetylation/Tests/test_", ti, ".results.tsv"),
            sep = "\t",
            col.names = TRUE
        )
    }

    # find all the TADs related to these tests
    tads_of_interest <- unlist(GRangesList(lapply(
        test_idx_of_interest,
        function(ti) {
            these_tads <- range(
                tads_gr[which(
                    apply(
                        tads,
                        1,
                        function(r) {
                            # only get the TADs related to this test index
                            ti %in% r$test_index
                        }
                    )
                )]
            )
            these_tads$test_index <- ti
            return(these_tads)
        }
    )))

    # find all the induced regions that map to these TADs
    regions_of_interest <- as.data.table(findOverlaps(regions_gr, tads_of_interest))
    regions_of_interest[, test_index := tads_of_interest[subjectHits]$test_index]
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
        by = test_index
    ]
    colnames(stouffer_z)[2] <- "z"
    # for any test that has all NAs, make sure it's listed as such (because the sum will be 0)
    tests_with_all_nas <- regions_of_interest[, all(is.na(p)), by = "test_index"][V1 == TRUE, test_index]
    stouffer_z[test_index %in% tests_with_all_nas, z := NA]
    # calculate the associated p-value
    stouffer_z[, p := 2 * pnorm(-abs(z))]
    # because we're only testing these regions, not the others, perform p-value adjustment here
    stouffer_z[, padj := p.adjust(p, method = "fdr", n = .N)]

    # save to table for later
    for (ti in test_idx_of_interest) {
        if (tests[test_index == ti, !is.na(z)]) {
            cat("This test has already been assigned: ", ti, "\n")
        }
        # variables for rows to get/set without performing the same calculation multiple times
        this_tests_ti = which(tests$test_index == ti)
        this_stouffer_ti = which(stouffer_z$test_index == ti)
        this_tad_ti = which(tads_of_interest$test_index == ti)
        tests[this_tests_ti]$z <- stouffer_z[this_stouffer_ti]$z
        tests[this_tests_ti]$p <- stouffer_z[this_stouffer_ti]$p
        tests[this_tests_ti]$padj <- stouffer_z[this_stouffer_ti]$padj
        tests[this_tests_ti]$chr <- as.character(seqnames(tads_of_interest)[this_tad_ti])
        tests[this_tests_ti]$start <- start(tads_of_interest[this_tad_ti]) - 1
        tests[this_tests_ti]$end <- end(tads_of_interest[this_tad_ti])
    }
}

# ==============================================================================
# Save data
# ==============================================================================
# convert lists back to comma-separated columns
tests$breakpoint_indices <- unlist(lapply(tests$breakpoint_indices, paste, collapse = ","))
tests$mutated_in <- unlist(lapply(tests$mutated_in, paste, collapse = ","))

# write the tables with columns in a particular order
fwrite(
    tests[, .SD, .SDcols = c(
        "chr", "start", "end",
        "test_index", "breakpoint_indices", "mutated_in",
        "n_mut", "n_nonmut",
        "z", "p", "padj"
    )],
    "sv-disruption-tests.acetylation.tsv",
    sep = "\t",
    col.names = TRUE
)
fwrite(sample_corrs, "acetylation-correlation.tsv", sep = "\t", col.names = TRUE)
