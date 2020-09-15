# ==============================================================================
# Meta
# ==============================================================================
# Calculate differential interactions
# --------------------------------------
# Description: Calculate differential interactions between benign and tumour samples
# Author: James Hawley


# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
suppressMessages(library("multiHiCcompare"))
suppressMessages(library("BiocParallel"))

# ==============================================================================
# Functions
# ==============================================================================
#' Read in contacts data from exported cooler dump file into format usable by multiHiCcompare
#'
#' @param path File path
#' @return sparse upper triangular matrix
read_contacts <- function(path) {
    dt <- fread(path, sep = "\t", header = FALSE)
    return(cooler2sparse(dt))
}

# ==============================================================================
# Data
# ==============================================================================
# read in metadata
metadata <- fread("config.tsv", sep = "\t", header = TRUE)
metadata <- metadata[Include == "Yes"]
SAMPLES <- metadata[, SampleID]
BENIGN_SAMPLES <- metadata[Source == "Primary" & Type == "Benign", SampleID]
TUMOUR_SAMPLES <- metadata[Source == "Primary" & Type == "Malignant", SampleID]

# load contacts
contacts <- lapply(
    SAMPLES,
    function(s) {
        read_contacts(file.path("TMP", paste0(s, ".300000000.res_40000bp.contacts.txt")))
    }
)
names(contacts) <- SAMPLES

# load blacklist and other regions to remove
data(hg38_cyto)

# ==============================================================================
# Analysis
# ==============================================================================
# register the number of cores to use simultaneously
num_cores <- 12
register(MulticoreParam(workers = num_cores), default = TRUE)

# create hicexp object
hic_expt <- make_hicexp(
    contacts,
    data_list = TRUE,
    groups = rep(c("Tumour", "Benign"), c(12, 5)),
    revome.regions = hg38_cyto
)

# free memory
rm(contacts)

# perform cyclical LOESS normalization
# (because we are using downsampled contacts, we don't need to adjust for library size)
hic_expt <- fastlo(hic_expt, verbose = TRUE, parallel = TRUE)

# plot fastlo normalization
png("Plots/fastlo-normalization.png", width = 12, height = 12, units = "cm", res = 300)
MD_hicexp(hic_expt)
dev.off()


# perform differential detection
hic_expt <- hic_exactTest(hic_expt, p.method = "fdr", parallel = TRUE)

# plot differntial MD plot
png("Plots/differential-md.png", width = 12, height = 12, units = "cm", res = 300)
MD_composite(hic_expt)
dev.off()

# extract differential results and save
res_dt <- data.table(results(hic_expt))
fwrite(res_dt, "difftl-interactions/difftl-results.tsv", sep = "\t", col.names = TRUE)

# free memory
rm(res_dt)

# make BED and BEDPE files of significantly differential interactions
loci <- data.table(topDirs(
    hic_expt,
    logfc_cutoff = 0.5,
    logcpm_cutoff = 0.5,
    p.adj_cutoff = 0.1,
    return_df = "bed"
))
fwrite(loci, "difftl-interactions/sig-loci.bed", sep = "\t", col.names = TRUE)


loci_pairs <- data.table(topDirs(
    hic_expt,
    logfc_cutoff = 0.5,
    logcpm_cutoff = 0.5,
    p.adj_cutoff = 0.1,
    return_df = "pairedbed"
))
fwrite(loci_pairs, "difftl-interactions/sig-loci-pairs.bedpe", sep = "\t", col.names = TRUE)
