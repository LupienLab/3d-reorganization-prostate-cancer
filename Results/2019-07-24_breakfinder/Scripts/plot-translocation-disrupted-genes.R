# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))

# ==============================================================================
# Data
# ==============================================================================
# load list of disrupted genes
genes <- fread(
    "translocation-disrupted-genes.tsv",
    sep = "\t",
    header = FALSE,
    col.names = c("chr", "start", "end", "name", "score", "strand", "EnsemblID")
)
# drop "score" column
genes[, score := NULL]


# load gene expression for 13 patients
exprs <- fread(
    file.path(
        "..", "..", "Data", "External", "CPC-GENE",
        "CPC-GENE_Chen-2019_RNAseq_rsem_gene_FPKM.13-Low-C-only.tsv"
    ),
    sep = "\t",
    header = TRUE
)

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

# calculate distance from insertion site
# new TAD boundary forms at chr14:35720000
new_boundary_pos <- 35720000
disrupted_exprs[
    chr == "chr14",
    distance_from := pmin(
        abs(start - new_boundary_pos),
        abs(end - new_boundary_pos)
    )
]

disrupted_exprs[, Non_Translocation_Mean := rowMeans(disrupted_exprs[, .SD, .SDcols = 2:13])]
disrupted_exprs[, Non_Translocation_SD := apply(
    disrupted_exprs[, .SD, .SDcols = 2:13],
    1,
    sd
)]
disrupted_exprs[, fold_change := PCa13848 / Non_Translocation_Mean]
disrupted_exprs[, z := (PCa13848 - Non_Translocation_Mean) / Non_Translocation_SD]

# perform t-test on disrupted_exprs$z
test_results <- t.test(disrupted_exprs$z)

# ==============================================================================
# Plots
# ==============================================================================
gg = (
    ggplot(data = disrupted_exprs_long)
    + geom_point(
        aes(x = (SampleID == "PCa13848"), y = Expression, colour = chr),
        position = "jitter"
    )
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
    height = 40,
    width = 24,
    units = "cm"
)

gg = (
    ggplot(data = disrupted_exprs)
    + geom_histogram(aes(x = z, y = ..density..))
    + stat_function(
        fun = dt,
        args = list(df = disrupted_exprs[, .N] - 1),
        linetype = 2
    )
    + stat_function(
        fun = dt,
        args = list(
            df = disrupted_exprs[, .N] - 1,
            ncp = test_results$estimate
        )
    )
    + geom_path(
        data = data.table(
            x = c(0, 0, test_results$estimate, test_results$estimate),
            y = c(0.4, 0.41, 0.41, 0.4)
        ),
        aes(x = x, y = y)
    )
    + annotate(
        geom = "text",
        x = test_results$estimate / 2,
        y = 0.41,
        vjust = -0.6,
        label = paste("p =", round(test_results$p.value, 4))
    )
    + labs(x = "z-score", y = "Density")
    + theme_minimal()
    + theme(
          axis.text.x = element_text(angle = 90)
      )
)
ggsave(
    "Plots/translocation-disrupted-genes.nhst.png",
    height = 12,
    width = 20,
    units = "cm"
)
