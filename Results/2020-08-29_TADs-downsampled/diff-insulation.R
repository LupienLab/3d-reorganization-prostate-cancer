# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("logging"))


# ==============================================================================
# Data
# ==============================================================================
loginfo("Loading Data")
# load metadata
metadata <- fread("config.tsv", sep = "\t", header = TRUE)
metadata <- metadata[Include == "Yes"]
SAMPLES <- metadata[, SampleID]
TUMOUR_SAMPLES <- metadata[Source == "Primary" & Type == "Malignant", SampleID]
BENIGN_SAMPLES <- metadata[Source == "Primary" & Type == "Benign", SampleID]
PRIMARY_SAMPLES <- metadata[Source == "Primary", SampleID]
n_tumour <- length(TUMOUR_SAMPLES)
n_benign <- length(BENIGN_SAMPLES)

# load TADs
insulation <- rbindlist(lapply(
    PRIMARY_SAMPLES,
    function(s) {
        dt <- fread(
            file.path("Insulation", paste0(s, ".300000000.insulation.tsv")),
            sep = "\t",
            header = TRUE
        )
        # remove "bad bins"
        dt <- dt[is_bad_bin == FALSE, .SD]
        dt[, SampleID := s]
        return(dt)
    }
))

# merge metadata with insulation scores
insulation <- merge(
    x = insulation,
    y = metadata[, .SD, .SDcols = c("SampleID", "Type")],
    by = "SampleID"
)


# ==============================================================================
# Analysis
# ==============================================================================
loginfo("Calculating differential insulation")
# calculate mean and SD of insulation at each region
diff_insulation <- insulation[,
    .(
        mean_ins = mean(2^log2_insulation_score_120000, na.rm = TRUE),
        sd_ins = sd(2^log2_insulation_score_120000, na.rm = TRUE)
    ),
    by = c("chrom", "start", "end", "Type")
]

# remove positions with no variance
diff_insulation <- diff_insulation[sd_ins > 0, .SD]
# re-cast table to have a single position for each row
wide_diff_ins <- dcast(diff_insulation,  chrom + start + end ~ Type, value.var = c("mean_ins", "sd_ins"))

# calculate pooled means and variance
wide_diff_ins[, `:=`(
    diff_mean_ins = mean_ins_Malignant - mean_ins_Benign,
    pooled_var = ((n_tumour - 1) * sd_ins_Malignant + (n_benign - 1) * sd_ins_Benign) / (n_tumour + n_benign - 2)
)]

# save results
fwrite(wide_diff_ins, "Insulation/difftl-insulation.tsv", sep = "\t", col.names = TRUE)

