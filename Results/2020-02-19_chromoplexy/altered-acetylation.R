# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("DESeq2"))
suppressMessages(library("GenomicRanges"))
suppressMessages(library("pheatmap"))
suppressMessages(library("RColorBrewer"))
suppressMessages(library("ggplot2"))

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
    splitv <- lapply(v, function(x) {strsplit(x, ",")[[1]]})
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

# load breakpoints
breakpoints <- fread("sv-breakpoints.tsv", sep = "\t", header = TRUE)

# load TADs encompassing each breakpoint
tads <- fread("sv-breakpoints.TADs.tsv", sep = "\t", header = TRUE)
tads_gr <- dt2gr(tads)

# load tests metadata
tests <- fread("sv-disruption-tests.tsv", sep = "\t", header = TRUE)
# convert comma-separated values into lists
tests$breakpoint_indices <- split_comma_col(tests$breakpoint_indices, as.numeric)
tests$mutated_in <- split_comma_col(tests$mutated_in)

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

# create design matrices for each set of mut-vs-nonmut (this depends on the breakpoint being considered)
all_comparisons <- unique(tests$mutated_in)

# perform differential analysis for each test
for (mut_samples in all_comparisons) {
    print(mut_samples)
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
    res <- as.data.table(results(dds))

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
    # find all the breakpoints related to these tests
    breakpoint_idx_of_interest <- tests[
        test_index %in% test_idx_of_interest,
        unique(unlist(breakpoint_indices))
    ]

    # find all the TADs related to these breakpoints
    tads_of_interest <- tads_gr[tads_gr$breakpoint_index %in% breakpoint_idx_of_interest]
    # find all the induced regions that map to these TADs
    regions_of_interest <- as.data.table(findOverlaps(regions_gr, tads_of_interest))
    # get differential analysis results for these regions
    regions_of_interest[, p := res[queryHits, pvalue]]
    regions_of_interest[, padj := res[queryHits, padj]]
    # calculate region weights for the TAD it is induced from
    regions_of_interest[, weight := sqrt(
        width(pintersect(
            regions_gr[queryHits],
            tads_of_interest[subjectHits]
        )) / width(tads_of_interest[subjectHits])
    )]
    # calculate Stouffer's Z with the weights calculated above
    stouffer_z <- regions_of_interest[
        !is.na(p) & !is.na(padj),
        sum(weight * qnorm(1 - p)) / sqrt(sum(weight ^ 2)),
        by = subjectHits
    ]
    colnames(stouffer_z)[2] <- "z"
    # calculate the associated p-value
    stouffer_z[, p := 2 * pnorm(-abs(z))]
    # because we're only testing these regions, not the others, perform p-value adjustment here
    stouffer_z[, padj := p.adjust(p, method = "fdr")]

    # map these values back to the `tads` object
    stouffer_z[, breakpoint_index := apply(
        .SD,
        1,
        function(r) tads_of_interest[r["subjectHits"]]$breakpoint_index
    )]
    tads <- merge(
        x = tads,
        y = stouffer_z[, .SD, .SDcols = -1],
        by = "breakpoint_index",
        all.x = TRUE
    )
}

# ==============================================================================
# Plots
# ==============================================================================
# heatmap of acetylation across the genome for all samples
ann_cols <- data.frame(
    T2E = metadata[, get("T2E Status")],
    LibrarySize = library_sizes
)
rownames(ann_cols) <- SAMPLES
pheatmap(
    mat = sample_corrs,
    color = colorRampPalette(brewer.pal(n = 7, name = "YlOrRd"))(100),
    #breaks = seq(0.8, 1, 0.01),
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    clustering_method = "ward.D2",
    annotation_col = ann_cols,
    legend = TRUE,
    filename = "Plots/H3K27ac-correlation-tad-induced.png"
)


# ==============================================================================
# Save data
# ==============================================================================
fwrite(tads, "sv-disruption-tests.acetylation.tsv", sep = "\t", col.names = TRUE)
