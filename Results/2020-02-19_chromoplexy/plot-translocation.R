# This script plots various stats about the T2E translocation in PCa13848

# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
source("plotting-helper.R")

GRAPH_DIR <- "Graphs"
# ==============================================================================
# Functions
# ==============================================================================

# ==============================================================================
# Data
# ==============================================================================
# load z-scores for tested genes
tested_genes <- fread(
    file.path(GRAPH_DIR, "sv-disruption-tests.expression.gene-level.tsv"),
    sep = "\t",
    header = TRUE
)

# load expression of all genes
exprs <- fread(
    file.path(
        "..",
        "..",
        "Data",
        "External",
        "CPC-GENE",
        "CPC-GENE_Chen-2019_RNAseq_rsem_gene_FPKM.13-LowC-only.tsv"
    ),
    sep = "\t",
    header = TRUE
)
# drop non-T2E, ETS+ sample
exprs[, PCa57054 := NULL]
exprs[, EnsemblID_short := gsub("\\.\\d+", "", EnsemblID)]

# load GENCODE reference annotation (all genes, not just protein-coding)
gencode <- fread(
    file.path("..", "..", "Data", "External", "GENCODE", "gencode.v33.all-genes.bed"),
    sep = "\t",
    header = FALSE,
    col.names = c("chr", "start", "end", "strand", "EnsemblID", "name"),
)
# remove annotation version number to ensure compatibility with previous RNA-seq
gencode[, EnsemblID_short := gsub("\\.\\d+", "", EnsemblID)]

# merge gencode annotations to gene expression
exprs <- merge(
    x = exprs,
    y = gencode,
    by = "EnsemblID_short"
)

# plotting thresholds
log_fold_thresh <- log2(c(0.5, 2))
abs_abundance_thresh <- c(-1, 1)

# delta offset used in calculating fold change
offset <- 1e-3

# chrom sizes
chrom_sizes <- fread(
    file.path("..", "..", "Data", "Processed", "2019-06-18_PCa-LowC-sequencing", "hg38.sizes.txt"),
    sep = "\t",
    header = FALSE,
    col.names = c("chrom", "size")
)

chrom_sizes[, colour := ifelse(.I %% 2 == 0, "#276FBF", "#183059")]
chrom_sizes[, offset := cumsum(as.numeric(c(0, head(size, -1))))]
chrom_sizes[, label_offset := (offset + c(0, head(offset, -1))) / 2]

# load sample metadata
metadata <- fread(
    file.path("..", "..", "Data", "External", "LowC_Samples_Data_Available.tsv"),
    sep = "\t",
    header = TRUE
)
metadata <- metadata[Include == "Yes"]
metadata[, SampleID := paste0("PCa", get("Sample ID"))]
SAMPLES <- metadata[, SampleID]
T2E_SAMPLES <- metadata[get("T2E Status") == "Yes", SampleID]
NONT2E_SAMPLES <- metadata[get("T2E Status") == "No", SampleID]
T2E_NONTRANSLOC_SAMPLES <- metadata[get("T2E Status") == "Yes" & SampleID != "PCa13848", SampleID]

# ==============================================================================
# Analysis
# ==============================================================================
# classify according to thresholds
tested_genes[, Fold_Thresh := (abs(log2fold) >= log_fold_thresh[2])]
tested_genes[, Abs_Thresh := (abs(mut_mean - nonmut_mean) >= abs_abundance_thresh[2])]
tested_genes[, Pass_Thresh := (Fold_Thresh & Abs_Thresh)]

# get genes associated with the translocation
insertion_site_genes <- tested_genes[test_ID == 44]
# transform into z-scores
insertion_site_genes_z <- insertion_site_genes[,
    .(chr, start, end, name, strand, EnsemblID_short, apply(.SD, 2, function(s) (s - nonmut_mean) / nonmut_sd)),
    .SDcols = colnames(insertion_site_genes)[grepl("PCa", colnames(insertion_site_genes))]
]
long_insertion_site_genes_z <- melt(
    insertion_site_genes_z,
    id.vars = c("chr", "start", "end", "name", "strand", "EnsemblID_short"),
    variable.name = "SampleID",
    value.name = "z"
)

# get genes in the translocated segment
ERG_END <- exprs[Symbol == "ERG", end]
TMPRSS2_START <- exprs[Symbol == "TMPRSS2", start]

translocated_genes <- exprs[chr == "chr21" & start >= ERG_END & end <= TMPRSS2_START]
translocated_genes[, log2fold := log2((PCa13848 + offset) / (apply(.SD, 1, mean) + offset)), .SDcols = T2E_NONTRANSLOC_SAMPLES]
translocated_genes[, mut_mean := PCa13848]
translocated_genes[, nonmut_mean := apply(.SD, 1, mean), .SDcols = T2E_NONTRANSLOC_SAMPLES]
translocated_genes[, Fold_Thresh := (abs(log2fold) >= log_fold_thresh[2])]
translocated_genes[, Abs_Thresh := (abs(mut_mean - nonmut_mean) >= abs_abundance_thresh[2])]
translocated_genes[, Pass_Thresh := (Fold_Thresh & Abs_Thresh)]

all_relevant_genes <- rbindlist(list(insertion_site_genes, translocated_genes), fill = TRUE)

# ==============================================================================
# Plots
# ==============================================================================
# plots associated with the T2E translocation from PCa13848
gg_insertion_site <- (
    ggplot()
    + geom_point(
        data = long_insertion_site_genes_z[SampleID != "PCa13848" & !is.na(z)],
        aes(
            # order columns according to z of mutated sample
            x = factor(
                name,
                ordered = TRUE,
                levels = long_insertion_site_genes_z[SampleID == "PCa13848"][order(z), unique(name)]
            ),
            y = z
        ),
        position = position_jitter(width = 0.2),
        colour = "#DBDBDB"
    )
    + geom_boxplot(
        data = long_insertion_site_genes_z[SampleID != "PCa13848" & !is.na(z)],
        aes(x = name, y = z),
        outlier.shape = NA,
        alpha = 0.7
    )
    + geom_point(
        data = long_insertion_site_genes_z[SampleID == "PCa13848" & !is.na(z)],
        aes(x = name, y = z, colour = abs(z) >= 1),
        colour = "#FF6347",
        size = 4,
        shape = "diamond"
    )
    + scale_y_continuous(
        limits = c(-2, 5.5),
        breaks = seq(-2, 5)
    )
    + labs(x = "Gene", y = expression(log[2] * " Expression Fold Change"))
    + coord_flip()
    + theme_minimal()
    + theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.minor = element_blank()
    )
)
savefig(gg_insertion_site, "Plots/sv-disruption/translocation.insertion-site", height = 12, width = 10)

gg_translocated_segment <- (
    ggplot(data = translocated_genes)
    + geom_point(aes(
        x = factor(start, ordered = TRUE, levels = sort(start)),
        y = log2fold,
        colour = Abs_Thresh
    ))
    + scale_x_discrete(
        breaks = translocated_genes[, start],
        labels = translocated_genes[, name]
    )
    + scale_y_continuous(
        limits = c(-10, 4),
        breaks = c(-10, -7, -4, -1, 0, 1, 4)
    )
    + scale_colour_manual(
        breaks = c(TRUE, FALSE),
        values = c("#bdbdbd", "#000000")
    )
    + labs(x = "Gene", y = expression(log[2] * " Expression Fold Change"))
    + theme_minimal()
    + theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.minor = element_blank()
    )
)
savefig(gg_translocated_segment, "Plots/sv-disruption/translocation.translocated-segment")

gg_detectable <- (
    ggplot(data = all_relevant_genes[Abs_Thresh == TRUE])
    + geom_histogram(aes(x = log2fold))
    # + geom_point(aes(x = 0, y = log2fold), position = position_jitter(height = 0, width = 0.2))
    # + geom_boxplot(aes(x = 0, y = log2fold), outlier.shape = NA, alpha = 0.2)
    + scale_x_continuous(
        limits = c(-4, 4)
    )
    + scale_colour_manual(
        breaks = c(TRUE, FALSE),
        values = c("#bdbdbd", "#000000")
    )
    + labs(y = "Gene Count", x = expression(log[2] * " Expression Fold Change"))
    + theme_minimal()
    + theme(
        panel.grid.minor = element_blank()
    )
)
savefig(gg_detectable, "Plots/sv-disruption/translocation.all-relevant-genes")
