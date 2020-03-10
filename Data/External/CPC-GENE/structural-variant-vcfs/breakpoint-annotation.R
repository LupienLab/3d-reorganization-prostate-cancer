# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
suppressMessages(library("GenomicRanges"))
suppressMessages(library("annotatr"))

# ==============================================================================
# Constants
# ==============================================================================

# ==============================================================================
# Data
# ==============================================================================
# load metadata
metadata <- fread("../../LowC_Samples_Data_Available.tsv", sep = "\t")
SAMPLES <- paste0("PCa", metadata[, get("Sample ID")])

# load SV breakpoints as GRanges objects
reg <- GRangesList(lapply(paste0(SAMPLES, ".hg38.bed"), read_regions, genome = "hg38"))
# add sample name as a metadata column
for (i in 1:length(reg)) {
    reg[[i]]$Sample = SAMPLES[i]
}

# ==============================================================================
# Analysis
# ==============================================================================
# annotations to extract
annots <- c(
    "hg38_genes_promoters",
    "hg38_genes_exons",
    "hg38_genes_introns",
    "hg38_genes_5UTRs",
    "hg38_genes_3UTRs",
    "hg38_genes_intergenic"
)

annotations <- build_annotations(genome = "hg38", annotations = annots)

# annotate the regions where SVs and their breakpoints are
reg_annotated <- GRangesList(lapply(
    reg,
    annotate_regions,
    annotations = annotations,
    ignore.strand = TRUE
))

# convert to data.table
regions_dt <- as.data.table(unlist(reg_annotated))

# reduce multiplicity of annotations, where possible
uniq_regions <- regions_dt[, .(unique(seqnames), unique(start) - 1, unique(end), unique(Sample), paste(unique(annot.type), sep = ",", collapse = ",")), by = c("name")]
colnames(uniq_regions) <- c("name", "chr", "start", "end", "sample", "annotation")

# save annotations in table
fwrite(uniq_regions, "sv-annotations.tsv", sep = "\t", col.names = TRUE)

# ==============================================================================
# Plots
# ==============================================================================