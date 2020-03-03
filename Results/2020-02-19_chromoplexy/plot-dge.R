# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))

# ==============================================================================
# Data
# ==============================================================================
# load hypothesis testing results
htest = fread("sv-disruption-tests.tsv", sep = "\t", header = TRUE)

# load z-scores for tested genes
tested_genes = fread("sv-disruption-tests.genes.tsv", sep = "\t", header = TRUE)

# load expression of all genes
exprs = fread(
    file.path(
        "..",
        "..",
        "Data",
        "External",
        "CPC-GENE",
        "CPC-GENE_Chen-2019_RNAseq_rsem_gene_FPKM.13-Low-C-only.tsv"
    ),
    sep = "\t",
    header = TRUE
)
exprs[, EnsemblID_short := gsub("\\.\\d+", "", EnsemblID)]

# load GENCODE reference annotation (all genes, not just protein-coding)
gencode = fread(
    file.path("..", "..", "Data", "External", "GENCODE", "gencode.v33.all-genes.bed"),
    sep = "\t",
    header = FALSE,
    col.names = c("chr", "start", "end", "strand", "EnsemblID", "name"),
)
# remove annotation version number to ensure compatibility with previous RNA-seq
gencode[, EnsemblID_short := gsub("\\.\\d+", "", EnsemblID)]

# merge gencode annotations to gene expression
exprs = merge(
    x = exprs,
    y = gencode,
    by = "EnsemblID_short"
)

# plotting thresholds
log_fold_thresh = log2(c(0.5, 2))
abs_abundance_thresh = c(-1, 1)

# delta offset used in calculating fold change
offset = 1e-3

# ==============================================================================
# Analysis
# ==============================================================================
# add n_breakpoints to tested_genes for plotting
tested_genes = merge(
    x = tested_genes,
    y = htest[, .SD, .SDcols = c("Sample", "Component_Index", "n_genes")],
    by.x = c("Mutated_In", "Component_Index"),
    by.y = c("Sample", "Component_Index"),
    all.x = TRUE
)

# calculate means and variances for all genes
sample_cols = 3:15
exprs[, Mean := rowMeans(.SD), .SDcols = sample_cols]
exprs[, StdDev := apply(.SD, 1, sd), .SDcols = sample_cols]

# calculate absolute abundances
tested_genes[, mutated_mean := (2^log2fold - 1) * (nonmut_mean + offset)]

# classify according to thresholds
tested_genes[, Pass_Thresh := FALSE]
tested_genes[abs(mutated_mean - nonmut_mean) >= abs_abundance_thresh[2] & abs(log2fold) >= log_fold_thresh[2], Pass_Thresh := TRUE]

# ==============================================================================
# Plots
# ==============================================================================
# p-value histogram across all connected components (i.e. chromoplexic events)
gg = (
    ggplot(data = htest)
    + geom_histogram(aes(x = p))
    + labs(x = "p-value", y = "Frequency")
    + theme_minimal()
)
ggsave(
    "Plots/sv-disruption.expression.p-values.png",
    height = 12,
    width = 20,
    units = "cm"
)

# plot z-scores for each connected component
gg = (
    ggplot(data = tested_genes)
    + geom_violin(aes(x = factor(Component_Index), y = z, fill = n_breakpoints))
    + labs(x = "Complex Event", y = "z-score\n(Mutated Sample)")
    + ylim(-5, 5)
    + facet_wrap(~ Mutated_In, scales = "free_y")
    + guides(fill = guide_legend(title="Breakpoints"))
    + scale_fill_viridis_c()
    + coord_flip()
    + theme_minimal()
    + theme(
        strip.text.y = element_text(angle = 0, hjust = 0),
        axis.text.y = element_blank(),
        legend.position = "bottom"
    )
)
ggsave(
    "Plots/sv-disruption.z.png",
    height = 20,
    width = 20,
    units = "cm"
)

gg = (
    ggplot(data = tested_genes)
    + geom_point(
        aes(x = factor(Component_Index), y = log2fold, group = factor(log2fold)),
        position = position_dodge(width=1),
        size = 1
    )
    # + labs(x = "Complex Event", y = expression(log[2] * " Expression fold change"))
    + ylim(tested_genes[, min(log2fold)], tested_genes[, -min(log2fold)])
    + facet_wrap(~ Mutated_In, scales = "free_y")
    + guides(fill = guide_legend(title="Breakpoints"))
    + coord_flip()
    + theme_minimal()
    + theme(
        strip.text.y = element_text(angle = 0, hjust = 0),
        axis.text.y = element_blank(),
        legend.position = "bottom"
    )
)
ggsave(
    "Plots/sv-disruption.fold-change.png",
    height = 30,
    width = 20,
    units = "cm"
)

gg = (
    ggplot(data = tested_genes)
    + geom_point(aes(x = mutated_mean, y = log2fold, colour = Pass_Thresh))
    + geom_hline(yintercept = log_fold_thresh[1], linetype = "dashed", colour = "red")
    + geom_hline(yintercept = log_fold_thresh[2], linetype = "dashed", colour = "red")
    + geom_vline(xintercept = abs_abundance_thresh[1], linetype = "dashed", colour = "red")
    + geom_vline(xintercept = abs_abundance_thresh[2], linetype = "dashed", colour = "red")
    # + labs(x = expression("RNA abundance difference ("*\bar{x}[mut] * ")"), y = expression(log[2] * " Expression fold change"))
    + labs(x = "RNA abundance difference", y = expression(log[2] * " Expression fold change"))
    + xlim(-100, 100)
    + scale_y_continuous(
        breaks = c(-20, -10, -5, -1, 0, 1, 5, 10)
    )
    + scale_colour_manual(
        limits = c(TRUE, FALSE),
        values = c("#000000", "#ececec")
    )
    + theme_minimal()
    + theme(
        legend.position = "bottom"
    )
)
ggsave(
    "Plots/sv-disruption.fold-change-vs-difference.png",
    height = 20,
    width = 20,
    units = "cm"
)

gg = (
    ggplot(data = exprs)
    + geom_density(aes(x = Mean))
    + labs(x = "Mean expression (FPKM)", y = "Density")
    # + xlim(0, 2000)
    + scale_x_log10()
    + theme_minimal()
    + theme(
        legend.position = "bottom"
    )
)
ggsave(
    "Plots/all-genes.sd-vs-mean.png",
    height = 20,
    width = 20,
    units = "cm"
)

# plot empirical CDFs of fold changes for each chromoxplexic event
# convert to log2 base for more understandable plotting
fold_ecdf = tested_genes[, ecdf(log2fold)]
ecdf_data = data.table(
    log2fold = tested_genes[order(log2fold), log2fold]
)
ecdf_data[, cdf := fold_ecdf(log2fold)]
gg = (
    ggplot()
    + geom_ribbon(
        data = ecdf_data[log2fold <= log_fold_thresh[1]],
        mapping = aes(x = log2fold, ymin = 0, ymax = cdf),
        alpha = 0.5,
        fill = "red"
    )
    + geom_ribbon(
        data = ecdf_data[log2fold >= log_fold_thresh[2]],
        mapping = aes(x = log2fold, ymin = cdf, ymax = 1),
        alpha = 0.5,
        fill = "red"
    )
    + stat_ecdf(
        data = tested_genes,
        mapping = aes(x = log2fold)
    )
    + geom_vline(xintercept=log_fold_thresh[1], linetype = "dashed", colour = "red")
    + geom_vline(xintercept=log_fold_thresh[2], linetype = "dashed", colour = "red")
    + annotate(
        x = -6, y = 0.25,
        label = paste0(round(fold_ecdf(log_fold_thresh[1]), 3) * 100, "%"),
        geom = "text",
        colour = "red",
        alpha = 0.5
    )
    + annotate(
        x = 4, y = 0.875,
        label = paste0((1 - round(fold_ecdf(log_fold_thresh[2]), 3)) * 100, "%"),
        geom = "text",
        colour = "red",
        alpha = 0.5
    )
    + geom_curve(
        aes(x = -5, y = 0.25, xend = -3, yend = 0.13),
        arrow = arrow(angle = 30, length = unit(0.03, "npc")),
        curvature = -0.3,
        colour = "red",
        alpha = 0.5
    )
    + geom_curve(
        aes(x = 3, y = 0.875, xend = 2, yend = 0.95),
        arrow = arrow(angle = 30, length = unit(0.03, "npc")),
        curvature = -0.3,
        colour = "red",
        alpha = 0.5
    )
    + scale_x_continuous(
        breaks = seq(-10, 10, 2),
        limits = c(-10, 10),
        name = expression(log[2] * " Expression Fold Change")
    )
    + scale_y_continuous(,
        breaks = seq(0, 1, 0.2),
        labels = seq(0, 100, 20),
        name = "Percentage of genes in TADs with SVs\n(Cumulative density)"
    )
    + theme_minimal()
)
ggsave(
    "Plots/sv-disruption.fold.ecdf.png",
    height = 12,
    width = 20,
    units = "cm"
)
