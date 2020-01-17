# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))

# ==============================================================================
# Data
# ==============================================================================
# load gencode data
gencode = fread(
    "gencode.v33.annotation.gtf",
    sep = "\t",
    skip = 5,
    header = FALSE,
    col.names = c("chr", "source", "feature", "start", "end", "score", "strand", "phase", "value")
)

# ==============================================================================
# Analysis
# ==============================================================================
# only keep protein coding genes
gencode_genes = gencode[feature == "gene" & grepl("protein_coding", value)]

# remove unnecessary columns
gencode_genes[, source := NULL]
gencode_genes[, phase := NULL]
gencode_genes[, feature := NULL]

# adjust start positions since they are 1-indexed
# see https://www.gencodegenes.org/pages/data_format.html
gencode_genes[, start := start - 1]

# create new columns from `value` column and parse out Ensembl IDs and gene names
value_data = gencode_genes[, tstrsplit(value, "; ", fixed = TRUE, keep = c(1, 3), names = c("Ensembl_ID", "Name"))]
value_data[, Ensembl_ID := gsub("gene_id ", "", Ensembl_ID)]
value_data[, Ensembl_ID := gsub("\"", "", Ensembl_ID)]
value_data[, Name := gsub("gene_name ", "", Name)]
value_data[, Name := gsub("\"", "", Name)]

gencode_genes[, value := NULL]
gencode_genes[, name := value_data$Name]
gencode_genes[, Ensembl_ID := value_data$Ensembl]

# ==============================================================================
# Save data
# ==============================================================================
# save in BED6 format
# see https://genome.ucsc.edu/FAQ/FAQformat.html#format1

fwrite(
    gencode_genes[, .SD, .SDcols= c("chr", "start", "end", "name", "score", "strand", "Ensembl_ID")],
    "gencode.v33.genes.bed",
    sep = "\t",
    col.names = FALSE
)
