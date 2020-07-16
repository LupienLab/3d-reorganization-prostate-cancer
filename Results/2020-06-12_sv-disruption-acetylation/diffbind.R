# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("DESeq2"))
suppressMessages(library("DiffBind"))
suppressMessages(library("GenomicRanges"))

BAM_DIR <- file.path("..", "..", "Data", "Processed", "2019-05-03_PCa-H3K27ac-peaks", "BAMs")
PEAK_DIR <- file.path("..", "..", "Data", "Processed", "2019-05-03_PCa-H3K27ac-peaks", "Peaks")
ACETYL_DIR <- file.path("Acetylation", "Peaks")
TEST_DIR <- file.path("Acetylation", "Tests")


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

# load TADs encompassing each breakpoint
tads <- fread(
    file.path("..", "2020-02-19_sv-disruption-TADs", "sv-disruption-tests.TADs.tsv"),
    sep = "\t",
    header = TRUE
)
tads_gr <- dt2gr(tads)

# load tests metadata
tests <- fread(
    file.path("..", "2020-02-19_chromoplexy", "Graphs", "sv-disruption-tests.tsv"),
    sep = "\t",
    header = TRUE
)
# convert comma-separated values into lists
tests$breakpoint_IDs <- split_comma_col(tests$breakpoint_IDs, as.numeric)
tests$mut_samples <- split_comma_col(tests$mut_samples)
tests$nonmut_samples <- split_comma_col(tests$nonmut_samples)

# load induced regions (the rows of the count matrix)
regions <- fread(
    file.path(PEAK_DIR, "catalogue-peaks.bed"),
    sep = "\t",
    header = FALSE,
    col.names = c("chr", "start", "end")
)
regions_gr <- dt2gr(regions)

# ==============================================================================
# Analysis
# ==============================================================================
# create design matrices for each set of mut-vs-nonmut (this depends on the test group)
# find all the unique combinations of foreground and background samples (mut vs nonmut)
# this will make creating the design matrices easier

all_comparisons <- unique(lapply(
    1:tests[, .N],
    function(i) {
        return(list(
            "mut" = tests[i, unlist(mut_samples)],
            "nonmut" = tests[i, unlist(nonmut_samples)]
        ))
    }
))

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

# remove all tests where a low quality sample is the mutated sample
all_comparisons_no_lowqual_mut <- Filter(
    Negate(is.null),
    lapply(
        all_comparisons,
        function(comp) {
            if ((length(comp$mut) == 1) && (comp$mut %in% low_qual_samples$SampleID)) {
                return(NULL)
            }
            return(comp)
        }
    )
)

# remove all low quality samples from remaining tests
all_comparisons_no_lowqual <- lapply(
    all_comparisons_no_lowqual_mut,
    function(comp) {
        return(list(
            mut = setdiff(comp$mut, low_qual_samples$SampleID),
            nonmut = setdiff(comp$nonmut, low_qual_samples$SampleID),
            test_IDs = comp$test_IDs
        ))
    }
)

# combine tests that, after removing low quality samples, now have the same split of samples
degenerate_comparisons <- unique(lapply(
    all_comparisons_no_lowqual,
    function(comp) {
        return(list(
            "mut" = comp$mut,
            "nonmut" = comp$nonmut
        ))
    }
))
reduced_comparisons_no_lowqual <- lapply(
    degenerate_comparisons,
    function(comp) {
        # get index of comparisons that have the same mut-vs-nonmut split
        comp_idx_mut <- which(sapply(all_comparisons_no_lowqual, function(all_comp) identical(all_comp$mut, comp$mut)))
        comp_idx_nonmut <- which(sapply(all_comparisons_no_lowqual, function(all_comp) identical(all_comp$nonmut, comp$nonmut)))
        comp_idx <- sort(unique(intersect(comp_idx_mut, comp_idx_nonmut)))
        combined_tests <- sort(unique(unlist(sapply(comp_idx, function(i) all_comparisons_no_lowqual[[i]]$test_IDs))))
        return(list(
            mut = comp$mut,
            nonmut = comp$nonmut,
            test_IDs = combined_tests
        ))
    }
)
all_comparisons <- reduced_comparisons_no_lowqual
rm(all_comparisons_no_lowqual_mut, all_comparisons_no_lowqual, degenerate_comparisons, reduced_comparisons_no_lowqual)


# list of regions of interest and the differential analysis results
roi <- list()

# load DBA object
dba_all <- readRDS("dba_all.rds")

# establish contrast for the specific comparison
for (i in 1:length(all_comparisons)) {
    cat("\n", i, "of", length(all_comparisons), "\n")
    comp <- all_comparisons[[i]]
    mut_samples <- comp$mut
    nonmut_samples <- comp$nonmut

    # establish contrast for differential analysis
    dba_comp <- dba.contrast(
        dba_all,
        group1 = metadata[, SampleID %in% mut_samples],
        group2 = metadata[, SampleID %in% nonmut_samples],
        name1 = "Mutated",
        name2 = "Nonmutated"
    )

    # perform differential analysis
    dba_comp <- dba.analyze(dba_comp)

    # extract differntially acetylated regions
    comp_regions <- dba.report(dba_comp, th = 1)
    # consider only the regions in the TADs around the breakpoints
    local_comp_regions <- subsetByOverlaps(
        comp_regions,
        tads[test_ID %in% comp$test_IDs, GRanges(
            seqnames = chr,
            ranges = IRanges(start = start, end = end)
        )]
    )
    # re-adjust q-values for only these peaks
    local_comp_regions$FDR <- p.adjust(local_comp_regions$`p-value`, method = "fdr")

    # convert results to data.table
    comp_dt <- as.data.table(comp_regions)
    comp_dt[, start := start - 1]
    colnames(comp_dt)[1] <- "chr"
    
    local_comp_dt <- as.data.table(local_comp_regions)
    local_comp_dt[, start := start - 1]
    colnames(local_comp_dt)[1] <- "chr"

    # save data for each relevant test
    for (ti in comp$test_IDs) {
        # write the full set of results
        fwrite(
            comp_dt,
            file.path(TEST_DIR, paste0("test_", ti, ".all.tsv")),
            sep = "\t",
            col.names = TRUE
        )

        # get loci only in the local region (i.e. the test TAD)
        test_tad <- tads[test_ID == ti]
        local_comp <- local_comp_dt[(chr == test_tad$chr) & (start <= test_tad$end) & (end >= test_tad$start)]

        fwrite(
            local_comp_dt,
            file.path(TEST_DIR, paste0("test_", ti, ".local.tsv")),
            sep = "\t",
            col.names = TRUE
        )
    }
}
