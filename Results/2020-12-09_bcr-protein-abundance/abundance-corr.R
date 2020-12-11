# ==============================================================================
# Meta
# ==============================================================================
# abundance-corr
# --------------------------------------
# Description: Test for the negative correlation between time to BCR and abundance of proteins involved in chromatin organization
# Author: James Hawley


# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("logging"))

loginfo("Loading packages")
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))

# ==============================================================================
# Functions
# ==============================================================================

# ==============================================================================
# Data
# ==============================================================================
loginfo("Loading data")
clinical <- fread("../../Data/External/CPC-GENE/CPC-GENE_Fraser-2017_Nature_Clinical-Data.tsv", sep = "\t", header = TRUE)
protein <- fread("../../Data/External/CPC-GENE/CPC-GENE_Sinha-Huang-2018_Proteomics_MassSpec.tsv", sep = "\t", header = TRUE)

structural_proteins <- c(
    "SMC1A",
    "SMC1B",
    "SMC3",
    "RAD21",
    "REC8",
    "STAG1",
    "STAG2",
    "STAG3",
    "PDS5A",
    "PDS5B",
    "WAPAL",
    "CDCA5",
    "ZNF143",
    "CTCF",
    "YY1"
)

protein <- protein[Gene %in% structural_proteins]

# load metadata for samples in our study
metadata <- fread("../../Data/External/LowC_Samples_Data_Available.tsv", sep = "\t", header = TRUE)
metadata <- metadata[Include == "Yes"]

# ==============================================================================
# Analysis
# ==============================================================================
loginfo("Performing calculations")

# convert to long form data table
protein_long <- melt(
    protein,
    id.vars = c("Protein Group", "Gene", "Protein name"),
    variable.name = "PatientID",
    value.name = "Abundance"
)
clinical[, PatientID := gsub("-F1", "", SampleID)]

# merge in time to BCR information
protein_long <- merge(
    x = protein_long,
    y = clinical[, .SD, .SDcols = c("PatientID", "BCR", "Time to BCR (years)")],
    by = "PatientID"
)

# test for correlation between "yes" and "no" BCR and each protein

protein_test <- protein_long[!is.na(BCR)]

cor_tests_BCR <- lapply(
    protein_test[, unique(Gene)],
    function(gene) {
        wilcox.test(
            x = protein_test[BCR == "No" & Gene == gene, Abundance],
            y = protein_test[BCR == "Yes" & Gene == gene, Abundance],
            alternative = "greater"
        )
    }
)
names(cor_tests_BCR) <- protein_test[, unique(Gene)]

# convert hypothesis test results into a data table
cor_tests_BCR_dt <- rbindlist(lapply(
    names(cor_tests_BCR),
    function(gene) {
        data.table(
            Gene = gene,
            W = cor_tests_BCR[[gene]]$statistic,
            p = cor_tests_BCR[[gene]]$p.value
        )
    }
))
cor_tests_BCR_dt[, FDR := p.adjust(p, method = "fdr")]

# test for time to BCR
cor_tests_BCR_time <- lapply(
    protein_test[, unique(Gene)],
    function(gene) {
        cor.test(
            x = protein_test[Gene == gene, Abundance],
            y = protein_test[Gene == gene, get("Time to BCR (years)")],
            alternative = "less",
            method = "spearman"
        )
    }
)
names(cor_tests_BCR_time) <- protein_test[, unique(Gene)]

# convert hypothesis test results into a data table
cor_tests_BCR_time_dt <- rbindlist(lapply(
    names(cor_tests_BCR_time),
    function(gene) {
        data.table(
            Gene = gene,
            rho = cor_tests_BCR_time[[gene]]$estimate,
            p = cor_tests_BCR_time[[gene]]$p.value
        )
    }
))
cor_tests_BCR_time_dt[, FDR := p.adjust(p, method = "fdr")]

# z-score transformation of abundances
gene_transform <- protein_long[,
    .(
        Mean = mean(Abundance, na.rm = TRUE),
        SD = sd(Abundance, na.rm = TRUE),
        Median = median(Abundance, na.rm = TRUE),
        IQR = quantile(Abundance, 0.75, na.rm = TRUE) - quantile(Abundance, 0.25, na.rm = TRUE)
    ),
    by = "Gene"
]

for (gene in protein_long[, unique(Gene)]) {
    mu <- gene_transform[Gene == gene, Mean]
    sigma <- gene_transform[Gene == gene, SD]
    protein_long[Gene == gene, Abundance_z := (Abundance - mu) / sigma]
}

protein_long[PatientID %in% metadata[, get("Patient ID")]]

# ==============================================================================
# Save data
# ==============================================================================
loginfo("Saving data")
fwrite(
    cor_tests_BCR_dt,
    "abundance-BCR-difference.tsv",
    sep = "\t",
    col.names = TRUE
)
fwrite(
    cor_tests_BCR_time_dt,
    "abundance-BCR-time-correlation.tsv",
    sep = "\t",
    col.names = TRUE
)
