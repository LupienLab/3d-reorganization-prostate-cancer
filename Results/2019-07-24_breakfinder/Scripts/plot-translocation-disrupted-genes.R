# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))

# ==============================================================================
# Data
# ==============================================================================
# load list of disrupted genes
genes = fread("translocation-disrupted-genes.tsv", sep = "\t", header = FALSE, col.names = c("chr", "start", "end", "name", "score", "strand", "EnsemblID"))
# drop "score" column
genes[, score := NULL]


# load gene expression for 13 patients
exprs = fread("../../Data/External/CPC-GENE/CPC-GENE_Chen-2019_RNAseq_rsem_gene_FPKM.13-Low-C-only.tsv", sep = "\t", header = TRUE)

# truncate Ensembl IDs to match across versions
genes[, EnsemblID := gsub("\\.\\d+$", "", EnsemblID)]
exprs[, EnsemblID := gsub("\\.\\d+$", "", EnsemblID)]

# combine the expression across patients for each gene in the dsirupted TAD
disrupted_exprs = merge(
    x = exprs,
    y = genes,
    by = "EnsemblID",
    all = FALSE
)

disrupted_exprs_long = melt(
    disrupted_exprs,
    id.vars = c("Symbol", "EnsemblID", "chr", "start", "end", "name", "strand"),
    variable.name = "SampleID",
    value.name = "Expression" 
)

# ==============================================================================
# Plots
# ==============================================================================
gg = (
    ggplot(data = disrupted_exprs_long)
    + geom_point(aes(x = (SampleID == "PCa13848"), y = Expression, colour = chr), position = "jitter")
    + labs(x = NULL, y = "Expression (FPKM)")
    + scale_x_discrete(
          limits = c(TRUE, FALSE),
          breaks = c(TRUE, FALSE),
          labels = c("PCa13848", "Other ERG+")
      )
    + facet_wrap(~ Symbol, scales = "free")
    + theme_minimal()
    + theme(
          axis.text.x = element_text(angle = 90)
      )
)
ggsave(
    "Plots/translocation-disrupted-genes.png",
    height = 20,
    width = 20,
    units = "cm"
)
