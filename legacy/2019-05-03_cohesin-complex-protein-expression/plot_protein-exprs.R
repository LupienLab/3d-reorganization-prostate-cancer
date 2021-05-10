# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))

# ==============================================================================
# Functions
# ==============================================================================
#' Save figures in multiple formats
#'
#' @param gg ggplot object
#' @param prefix Prefix for output file
#' @param ext Output extensions
#' @param dpi DPI resolution
savefig = function(gg, prefix, ext = c("png", "pdf"), width = 20, height = 12, dpi = 400) {
    for (e in ext) {
        ggsave(
            paste(prefix, e, sep = "."),
            gg,
            height = height,
            width = width,
            units = "cm",
            dpi = dpi
        )
    }
}

# ==============================================================================
# Data
# ==============================================================================
# metadata of patient samples
metadata <- fread(
    "../../Data/External/LowC_Samples_Data_Available.tsv",
    sep = "\t",
    header = TRUE
)
# only keep included samples
metadata <- metadata[Include == "Yes"]

# change header for simpler acces
colnames(metadata) = gsub("[ -]", "_", colnames(metadata))

# protein information from Sinha, Huang et al
protein_exprs = fread(
    "../../Data/External/CPC-GENE/CPC-GENE_Sinha-Huang-2018_Proteomics_MassSpec.tsv",
    sep = "\t",
    header = TRUE
)
# subset data to only consider patients in our Low-C dataset
col_idx = c(1:3, which(colnames(protein_exprs) %in% metadata$Patient_ID))
protein_exprs = protein_exprs[, ..col_idx]

# get list of proteins to compare expression for
test_proteins = fread("test-proteins.tsv", header = TRUE)
# only keep proteins that exist in protein expression data
test_proteins = test_proteins[which(Name %in% protein_exprs$Gene), ]

patient_exprs = melt(
    protein_exprs[Gene %in% test_proteins$Name, .SD],
    id.vars = c("Protein Group", "Gene", "Protein name"),
    variable.name = "Patient_ID",
    value.name = "Protein_Level"
)

# ==============================================================================
# Analysis
# ==============================================================================
# calculate the median absolute deviaton for all proteins in these 13 patients
mat = as.matrix(protein_exprs[, 4:dim(protein_exprs)[2]])
plotting_stats = data.table(
    Gene = protein_exprs[, Gene],
    Median = apply(mat, 1, median, na.rm = TRUE),
    MAD = apply(mat, 1, mad, na.rm = TRUE)
)
q95_mad_all_proteins = quantile(plotting_stats$MAD, 0.95, na.rm = TRUE)

# add plotting stats to patient_exprs
patient_exprs = merge(
    patient_exprs,
    plotting_stats,
    by = "Gene"
)
patient_exprs[, Lower := Median - q95_mad_all_proteins]
patient_exprs[, Upper := Median + q95_mad_all_proteins]

# ==============================================================================
# Plots
# ==============================================================================
gg <- (
    ggplot(data = patient_exprs, mapping = aes(x = Patient_ID))
    + geom_col(aes(y = Protein_Level, fill = Patient_ID))
    # + geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = Lower, ymax = Upper), alpha = 0.1)
    # + geom_hline(aes(yintercept =  Median), linetype = "dashed")
    + labs(x = NULL, y = expression(log[2] * "(Protein Intensity)"))
    + scale_x_discrete(
        limits = metadata[order(Patient_ID), Patient_ID],
        labels = metadata[order(Patient_ID), Patient_ID]
    )
    + scale_fill_manual(
        limits = metadata[, Patient_ID],
        values = metadata[, Colour]
    )
    + facet_wrap(~ Gene)
    + guides(fill = FALSE)
    + theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    )
)
savefig(gg, "loop-forming-factors_expression")
