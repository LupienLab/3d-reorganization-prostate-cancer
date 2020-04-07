# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))

GRAPH_DIR <- "Graphs"

# ==============================================================================
# Data
# ==============================================================================
# load hypothesis testing results
htest <- fread(file.path(GRAPH_DIR, "sv-disruption-tests.expression.tsv"), sep = "\t", header = TRUE)

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
        "CPC-GENE_Chen-2019_RNAseq_rsem_gene_FPKM.13-Low-C-only.tsv"
    ),
    sep = "\t",
    header = TRUE
)
exprs[, EnsemblID_short := gsub("\\.\\d+", "", EnsemblID)]

# load GENCODE reference annotation (all genes, not just protein-coding)
gencode <- fread(
    file.path(
        "..", "..", "Data", "External",
        "GENCODE", "gencode.v33.all-genes.bed"
    ),
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
    "../../Data/Processed/2019-06-18_PCa-LowC-sequencing/hg38.sizes.txt",
    sep = "\t",
    header = FALSE,
    col.names = c("chrom", "size")
)

chrom_sizes[, colour := ifelse(.I %% 2 == 0, "#276FBF", "#183059")]
chrom_sizes[, offset := cumsum(as.numeric(c(0, head(size, -1))))]
chrom_sizes[, label_offset := (offset + c(0, head(offset, -1))) / 2]

# ==============================================================================
# Analysis
# ==============================================================================
# get IDs of test groups where there is significantly altered gene expression
altering_tests <- htest[FDR < 0.05, test_ID]

# calculate means and variances for all genes
sample_cols <- 3:15
exprs[, means := rowMeans(.SD), .SDcols = sample_cols]
exprs[, std_dev := apply(.SD, 1, sd), .SDcols = sample_cols]

# classify according to thresholds
tested_genes[, Pass_Thresh := FALSE]
tested_genes[
    (
        abs(mut_mean - nonmut_mean) >= abs_abundance_thresh[2]
        & abs(log2fold) >= log_fold_thresh[2]
    ),
    Pass_Thresh := TRUE
]

# order chromosomes
CHRS <- paste0("chr", c(1:22, "X", "Y"))
tested_genes[, chr := factor(chr, ordered = TRUE, levels = CHRS)]

# set cumulative genome position for Manhattan-style plotting
tested_genes$cpos <- sapply(
    1:tested_genes[, .N],
    function(i) {
        return(
            tested_genes[i, start]
            + chrom_sizes[chrom == tested_genes[i, chr], offset]
        )
    }
)
tested_genes$colour <- sapply(
    1:tested_genes[, .N],
    function(i) {
        return(chrom_sizes[chrom == tested_genes[i, chr], colour])
    }
)

# state whether the genes are related to a breakpoint that significantly
# alters the set of genes (doing this calculation per row with ifelse)
tested_genes[, altering_test := ifelse(
    test_ID %in% altering_tests,
    TRUE,
    FALSE
)]

# set of tests where at least 1 gene passes both thresholds
tests_with_some_deg <- tested_genes[
    Pass_Thresh == TRUE,
    .N,
    by = test_ID
]$test_ID

# calculate percentage of genes in these TADs
tests_altered_genes <- data.table(
    test_ID = tests_with_some_deg,
    n_altered_genes = tested_genes[
        Pass_Thresh == TRUE,
        .N,
        by = test_ID
    ]$N,
    n_genes = tested_genes[
        test_ID %in% tests_with_some_deg,
        .N,
        by = test_ID
    ]$N
)
tests_altered_genes[, pct_genes := n_altered_genes / n_genes]
print(tests_altered_genes[, summary(pct_genes)])

# ==============================================================================
# Plots
# ==============================================================================
# p-value histogram across all connected components (i.e. chromoplexic events)
gg_pval_hist <- (
    ggplot(data = htest)
    + geom_histogram(aes(x = p))
    + labs(x = "p-value", y = "Frequency")
    + theme_minimal()
)
ggsave(
    "Plots/sv-disruption/expression.p-values.png",
    gg_pval_hist,
    height = 12,
    width = 20,
    units = "cm"
)
ggsave(
    "Plots/sv-disruption/expression.p-values.pdf",
    gg_pval_hist,
    height = 12,
    width = 20,
    units = "cm"
)

# plot z-scores for each breakpoint
ann_text <- tested_genes[z > 20 & is.finite(z)]

gg_z <- (
    ggplot(data = tested_genes[is.finite(z)])
    + geom_text(
        data = ann_text,
        aes(x = cpos, y = z + 20, label = name),
    )
    + geom_point(
        aes(x = cpos, y = z, colour = colour),
        size = 1
    )
    + labs(x = "Breakpoint Position", y = "z-score")
    + scale_x_continuous(
        limits = c(0, chrom_sizes[chrom == "chrM", offset]),
        breaks = chrom_sizes[, label_offset],
        label = chrom_sizes[, chrom]
    )
    + guides(colour = FALSE)
    + theme_minimal()
    + theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 0.5),
        strip.text.x = element_text(angle = 90, hjust = 0.5),
        legend.position = "bottom"
    )
)
ggsave(
    "Plots/sv-disruption/expression.z.png",
    gg_z,
    height = 10,
    width = 40,
    units = "cm"
)
ggsave(
    "Plots/sv-disruption/expression.z.pdf",
    gg_z,
    height = 10,
    width = 40,
    units = "cm"
)

gg_fc <- (
    ggplot(data = tested_genes)
    + geom_boxplot(
        aes(x = factor(test_ID), y = log2fold, colour = altering_test)
    )
    + labs(x = "Test", y = expression(log[2] * " Expression fold change"))
    + ylim(tested_genes[, min(log2fold)], tested_genes[, -min(log2fold)])
    + guides(fill = guide_legend(title = "Related breakpoints"))
    + scale_colour_manual(
        limits = c(TRUE, FALSE),
        values = c("#000000", "#cecece"),
        name = "Breakpoints significantly alters expression"
    )
    + coord_flip()
    + theme_minimal()
    + theme(
        strip.text.y = element_text(angle = 0, hjust = 0),
        axis.text.y = element_blank(),
        legend.position = "bottom"
    )
)
ggsave(
    "Plots/sv-disruption/expression.fold-change.png",
    gg_fc,
    height = 30,
    width = 20,
    units = "cm"
)
ggsave(
    "Plots/sv-disruption/expression.fold-change.pdf",
    gg_fc,
    height = 30,
    width = 20,
    units = "cm"
)

# same as above, but thresholding on absolute difference in mRNA abundance
gg_fc_thresh <- (
    ggplot(data = tested_genes[
        abs(mut_mean - nonmut_mean) >= abs_abundance_thresh[2],
        .SD
    ])
    + geom_point(
        aes(
            x = factor(test_ID),
            y = log2fold,
            group = factor(log2fold),
            colour = (test_ID %% 2 == 0)
        ),
        position = position_dodge(width = 1),
        size = 1
    )
    + labs(x = "Test", y = expression(log[2] * " Expression fold change"))
    + ylim(tested_genes[, min(log2fold)], tested_genes[, -min(log2fold)])
    + guides(fill = FALSE)
    + coord_flip()
    + theme_minimal()
    + theme(
        strip.text.y = element_text(angle = 0, hjust = 0),
        axis.text.y = element_blank(),
        legend.position = "bottom"
    )
)
ggsave(
    "Plots/sv-disruption/expression.fold-change.thresholded.png",
    gg_fc_thresh,
    height = 30,
    width = 20,
    units = "cm"
)
ggsave(
    "Plots/sv-disruption/expression.fold-change.thresholded.pdf",
    gg_fc_thresh,
    height = 30,
    width = 20,
    units = "cm"
)

gg_fc_diff <- (
    ggplot(data = tested_genes)
    + geom_point(aes(
        x = mut_mean - nonmut_mean,
        y = log2fold,
        colour = Pass_Thresh
    ))
    + geom_hline(
        yintercept = log_fold_thresh[1],
        linetype = "dashed",
        colour = "red"
    )
    + geom_hline(
        yintercept = log_fold_thresh[2],
        linetype = "dashed",
        colour = "red"
    )
    + geom_vline(
        xintercept = abs_abundance_thresh[1],
        linetype = "dashed",
        colour = "red"
    )
    + geom_vline(
        xintercept = abs_abundance_thresh[2],
        linetype = "dashed",
        colour = "red"
    )
    + labs(
        x = "RNA abundance difference",
        y = expression(log[2] * " Expression fold change")
    )
    + xlim(-100, 100)
    + scale_y_continuous(
        breaks = c(-20, -10, -5, -1, 0, 1, 5, 10)
    )
    + scale_colour_manual(
        limits = c(TRUE, FALSE),
        values = c("#000000", "#ececec")
    )
    + guides(colour = FALSE)
    + theme_minimal()
    + theme(
        legend.position = "bottom"
    )
)
ggsave(
    "Plots/sv-disruption/expression.fold-change-vs-difference.png",
    gg_fc_diff,
    height = 20,
    width = 20,
    units = "cm"
)
ggsave(
    "Plots/sv-disruption/expression.fold-change-vs-difference.pdf",
    gg_fc_diff,
    height = 20,
    width = 20,
    units = "cm"
)

gg_genes <- (
    ggplot(data = exprs)
    + geom_density(aes(x = means))
    + labs(x = "Mean expression (FPKM)", y = "Density")
    + scale_x_log10()
    + theme_minimal()
    + theme(
        legend.position = "bottom"
    )
)
ggsave(
    "Plots/all-genes.sd-vs-mean.png",
    gg_genes,
    height = 20,
    width = 20,
    units = "cm"
)
ggsave(
    "Plots/all-genes.sd-vs-mean.pdf",
    gg_genes,
    height = 20,
    width = 20,
    units = "cm"
)

# plot empirical CDFs of fold changes for each breakpoint
fold_ecdf <- tested_genes[, ecdf(log2fold)]
ecdf_data <- data.table(
    log2fold = tested_genes[order(log2fold), log2fold]
)
ecdf_data[, cdf := fold_ecdf(log2fold)]
gg_fc_cdf <- (
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
    + geom_path(
        data = ecdf_data,
        aes(x = log2fold, y = cdf)
    )
    + geom_vline(
        xintercept = log_fold_thresh[1],
        linetype = "dashed",
        colour = "red"
    )
    + geom_vline(
        xintercept = log_fold_thresh[2],
        linetype = "dashed",
        colour = "red"
    )
    + annotate(
        x = -6, y = 0.25,
        label = paste0(round(fold_ecdf(log_fold_thresh[1]), 3) * 100, "%"),
        geom = "text",
        colour = "red",
        alpha = 0.5
    )
    + annotate(
        x = 4, y = 0.875,
        label = paste0(
            (1 - round(fold_ecdf(log_fold_thresh[2]), 3)) * 100, "%"
        ),
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
    "Plots/sv-disruption/expression.fold-change.ecdf.png",
    gg_fc_cdf,
    height = 12,
    width = 20,
    units = "cm"
)
ggsave(
    "Plots/sv-disruption/expression.fold-change.ecdf.pdf",
    gg_fc_cdf,
    height = 12,
    width = 20,
    units = "cm"
)

fold_ecdf <- tested_genes[
    abs(mut_mean - nonmut_mean) >= abs_abundance_thresh[2],
    ecdf(log2fold)
]
ecdf_data <- data.table(
    log2fold = tested_genes[
        abs(mut_mean - nonmut_mean) >= abs_abundance_thresh[2]
    ][order(log2fold), log2fold]
)
ecdf_data[, cdf := fold_ecdf(log2fold)]
gg_fc_thresh_cdf <- (
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
    + geom_path(
        data = ecdf_data,
        aes(x = log2fold, y = cdf)
    )
    + geom_vline(
        xintercept = log_fold_thresh[1],
        linetype = "dashed",
        colour = "red"
    )
    + geom_vline(
        xintercept = log_fold_thresh[2],
        linetype = "dashed",
        colour = "red"
    )
    + annotate(
        x = -6, y = 0.25,
        label = paste0(round(fold_ecdf(log_fold_thresh[1]), 3) * 100, "%"),
        geom = "text",
        colour = "red",
        alpha = 0.5
    )
    + annotate(
        x = 4, y = 0.875,
        label = paste0(
            (1 - round(fold_ecdf(log_fold_thresh[2]), 3)) * 100,
            "%"
        ),
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
    "Plots/sv-disruption/expression.fold-change.ecdf.thresholded.png",
    gg_fc_thresh_cdf,
    height = 12,
    width = 20,
    units = "cm"
)
ggsave(
    "Plots/sv-disruption/expression.fold-change.ecdf.thresholded.pdf",
    gg_fc_thresh_cdf,
    height = 12,
    width = 20,
    units = "cm"
)
