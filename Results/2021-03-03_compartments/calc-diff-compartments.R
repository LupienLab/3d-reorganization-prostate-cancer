# ==============================================================================
# Meta
# ==============================================================================
# compartments
# ------------------------------------------------
# Author: James Hawley
# Description: Compartments analysis


# ==============================================================================
# Environment
# ==============================================================================
suppressWarnings(library("logging"))
loginfo("Loading packages")

suppressWarnings(library("data.table"))
suppressWarnings(library("ggplot2"))

CMPMT_DIR <- file.path(
    "..", "..", "Data", "Processed", "2019-06-18_PCa-LowC-sequencing", "Contacts"
)

# ==============================================================================
# Functions
# ==============================================================================
lm2dt <- function(mod) {
    mod_summ <- summary(mod)
    res <- data.table(
        mu = mod_summ$coefficients["(Intercept)", "Estimate"],
        beta = mod_summ$coefficients["TypeMalignant", "Estimate"],
        beta_se = mod_summ$coefficients["TypeMalignant", "Std. Error"],
        t = mod_summ$coefficients["TypeMalignant", "t value"],
        df = mod$df.residual,
        p = mod_summ$coefficients["TypeMalignant", "Pr(>|t|)"]
    )
    return(res)
}


# ==============================================================================
# Data
# ==============================================================================
loginfo("Loading data")

# load sample metadata
meta <- fread("config.tsv")
meta <- meta[Include == "Yes"]
SAMPLES <- meta[, Sample_ID]

# load compartment eigenvalues
eigs <- rbindlist(lapply(
    SAMPLES,
    function(s) {
        dt <- fread(
            file.path(CMPMT_DIR, paste0(s, ".compartments.cis.vecs.tsv")),
            sep = "\t",
            header = TRUE
        )
        dt[, bin_ID := .I]
        dt[, Sample_ID := s]
        return(dt[,
            .SD,
            .SDcols = c("Sample_ID", "chrom", "start", "end", "bin_ID", "E1", "E2", "E3")
        ])
    }
))

# add tissue type metadata to eigs
eigs <- merge(
    x = eigs,
    y = meta[, .SD, .SDcols = c("Sample_ID", "Label", "Type")],
    by = "Sample_ID"
)

eigs[, Sample_ID := as.factor(Sample_ID)]
eigs[, bin_ID := as.factor(bin_ID)]
eigs[, Type := as.factor(Type)]
eigs <- eigs[complete.cases(eigs)]

# ==============================================================================
# Analysis
# ==============================================================================
loginfo("Calculating basic statistics")

# get unique bins and positions for future mapping
# don't need them for most of the calculations, so using a bin ID is more memory efficient
genomic_bins <- unique(eigs[, .SD, .SDcols = c("chrom", "start", "end", "bin_ID")])

# calculate mean, sd, and n of E1 for each bin
long_e1 <- eigs[,
    .(
        Mean = mean(E1, na.rm = TRUE),
        SD = sd(E1, na.rm = TRUE),
        N = .N
    ),
    by = c("bin_ID", "Type")
]

# cast to wide format for each variable
wide_e1_n <- dcast(
    long_e1,
    bin_ID ~ Type,
    value.var = "N"
)
wide_e1_mean <- dcast(
    long_e1,
    bin_ID ~ Type,
    value.var = "Mean"
)
wide_e1_sd <- dcast(
    long_e1,
    bin_ID ~ Type,
    value.var = "SD"
)

# combine each wide table into a single wide one
wide_e1 <- merge(
    x = wide_e1_n,
    y = wide_e1_mean,
    by = "bin_ID",
    suffixes = c("_N", "_Mean")
)
wide_e1 <- merge(
    x = wide_e1,
    # gracefully rename columns before merging
    y = wide_e1_sd[, .(bin_ID, Benign_SD = Benign, Malignant_SD = Malignant)],
    by = "bin_ID"
)
wide_e1 <- merge(
    x = wide_e1,
    y = eigs[,
        .(
            # gracefully rename columns before merging
            All_Mean = mean(E1, na.rm = TRUE),
            All_SD = sd(E1, na.rm = TRUE),
            All_N = .N
        ),
        by = "bin_ID"
    ],
    by = "bin_ID"
)

loginfo("Calculating linear models")
# check which bins have enough samples to calculate a contrast
design <- dcast(
    eigs[, .N, by = c("Type", "bin_ID")],
    bin_ID ~ Type,
    value.var = "N",
    fill = 0
)
# only keep bins with values from >= 1 benign and >= 1 malignant samples
design <- design[(Benign > 0) & (Malignant > 0), .SD]
max_bin_ID <- design[, max(as.integer(bin_ID))]

# calculating linear models for compartment differences
model <- rbindlist(lapply(
    design[, bin_ID],
    function(b) {
        cat(paste(as.integer(b), "of", max_bin_ID, "\r"))
        mod <- lm(E1 ~ Type, data = eigs[bin_ID == b])
        dt <- lm2dt(mod)
        dt[, bin_ID := b]
        return(dt)
    }
))
cat("\n")
model[, q := p.adjust(p, method = "fdr")]

# map back to genomic coordinates before saving
model <- merge(
    x = genomic_bins,
    y = model,
    by = "bin_ID",
    # keep all bins, just enforce NAs if the bin has been filtered out for any reason
    all = TRUE
)
wide_e1 <- merge(
    x = genomic_bins,
    y = wide_e1,
    by = "bin_ID",
    all = TRUE
)

wide_e1[, `:=`(
    Diff = Malignant_Mean - Benign_Mean,
    Abs_Diff = abs(Malignant_Mean - Benign_Mean)
)]

s <- wide_e1[!is.na(Abs_Diff) & (All_N > 15), .SD, keyby = -Abs_Diff]
s_sorted <- wide_e1[
    !(chrom %in% c("chr3", "chr8", "chr13", "chr19", "chrY")) & abs(Abs_Diff) > 0.5,
    .SD,
    keyby = "bin_ID"
]

# ==============================================================================
# Save Data
# ==============================================================================
loginfo("Saving data")
# save model outputs without the bin_ID
fwrite(
    model[, .SD, .SDcols = -"bin_ID"],
    "compartments.results.tsv",
    sep = "\t",
    col.names = TRUE
)

# save wide compartment data for other plotting purposes
fwrite(
    wide_e1[, .SD, .SDcols = -"bin_ID"],
    "compartments.stats.tsv",
    sep = "\t",
    col.names = TRUE
)
