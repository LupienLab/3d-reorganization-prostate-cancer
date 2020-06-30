# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
suppressMessages(library("sleuth"))
suppressMessages(library("gridExtra"))
source("../2020-02-19_chromoplexy/plotting-helper.R")

QVAL_THRESH <- 0.05
LOG2FOLD_THRESH <- 1

PLOT_DIR <- "Plots"


# ==============================================================================
# Functions
# ==============================================================================
split_comma_col <- function(v, f=identity) {
    # split into list
    splitv <- lapply(v, function(x) {strsplit(x, ",")[[1]]})
    # remove various non-informative characters (spaces, braces)
    splitv <- lapply(splitv, function(x) {gsub("[][ ]", "", x)})
    return(lapply(splitv, f))
}


# ==============================================================================
# Data
# ==============================================================================
# load tests
sv_tests <- fread(
    "../2020-02-19_chromoplexy/Graphs/sv-disruption-tests.tsv",
    sep = "\t",
    header = TRUE
)
sv_tests$mut_samples <- split_comma_col(sv_tests$mut_samples)
sv_tests$nonmut_samples <- split_comma_col(sv_tests$nonmut_samples)

SAMPLES <- sv_tests[, sort(unique(unlist(mut_samples)))]

# load TAD disruption test results
tad_tests <- fread(
    file.path("..", "2020-02-19_sv-disruption-TADs", "sv-disruption-tests.TADs.tsv"),
    sep = "\t",
    header = TRUE
)

gencode <- fread(
    file.path("..", "..", "Data", "External", "GENCODE", "gencode.v33.all-transcripts.bed"),
    sep = "\t",
    header = FALSE,
    col.names = c("chr", "start", "end", "strand", "gene_id", "gene_name", "target_id", "transcript_name")
)

# load all naive test results
tested_genes <- rbindlist(lapply(
    sv_tests$test_ID,
    function(tid) {
        dt <- fread(
            paste0("sleuth/test_", tid, ".genes.tested.tsv"),
            sep = "\t",
            header = TRUE
        )
        dt[, test_ID := tid]
        return(dt)
    }
))

tested_transcripts <- rbindlist(lapply(
    sv_tests$test_ID,
    function(tid) {
        dt <- fread(
            paste0("sleuth/test_", tid, ".transcripts.tested.tsv"),
            sep = "\t",
            header = TRUE
        )
        dt[, test_ID := tid]
        return(dt)
    }
))

transcript_table <- rbindlist(lapply(
    SAMPLES,
    function(s) {
        dt <- fread(
            paste0("../../Data/Processed/2020-06-17_PCa-RNA-seq/", s, "/abundance.tsv"),
            sep = "\t",
            header = TRUE
        )
        dt[, SampleID := s]
        return(dt)
    }
))


# ==============================================================================
# Analysis
# ==============================================================================
# combine transcripts for gene counts
gene_tpm <- rbindlist(apply(
    tested_genes,
    1,
    function(r) {
        # get gene_id
        gene <- r["target_id"]
        # get all transcripts associated with that gene
        transcript_ids <- gencode[gene_id == gene, target_id]
        # sum all TPM for all transcripts for that gene
        gene_tpm <- transcript_table[target_id %in% transcript_ids, .(tpm = sum(tpm)), by = "SampleID"]
        gene_tpm[, gene_id := gene]
        
        # get test_ID
        tid <- as.integer(r["test_ID"])
        # add Mutated column for each patient involved in a test
        gene_tpm[, test_ID := tid]
        gene_tpm[, Mutated := ""]
        mut_samples <- sv_tests[test_ID == tid, unlist(mut_samples)]
        nonmut_samples <- sv_tests[test_ID == tid, unlist(nonmut_samples)]
        gene_tpm[SampleID %in% mut_samples, Mutated := "Mutated"]
        gene_tpm[SampleID %in% nonmut_samples, Mutated := "Nonmutated"]

        # remove NA samples that are excluded from each test
        gene_tpm <- gene_tpm[Mutated != ""]
        return(gene_tpm)
    }
))

gene_tpm_mutated <- gene_tpm[, .(tpm = mean(tpm)), by = c("test_ID", "gene_id", "Mutated")]
gene_tpm_mutated_per_test <- dcast(gene_tpm_mutated, test_ID + gene_id ~ Mutated)
gene_tpm_mutated_per_test[, log2fc := log2((Mutated + 1) / (Nonmutated + 1))]

# merge with p-/q-values
merged_gene_tpm <- merge(
    x = tested_genes,
    y = gene_tpm_mutated_per_test,
    by.x = c("test_ID", "target_id"),
    by.y = c("test_ID", "gene_id")
)
merged_gene_tpm[, altered_group := (qval < QVAL_THRESH & abs(log2fc) > log(2))]
merged_gene_tpm[, test_ID := factor(test_ID)]

counted_gene_tpm <- rbindlist(lapply(
    sv_tests$test_ID,
    function(tid) {
        dt <- data.table(
            test_ID = tid,
            Significant = c(FALSE, TRUE, TRUE),
            LargeFold = c(NA, FALSE, TRUE),
            N = c(
                merged_gene_tpm[(test_ID == tid) & (qval > QVAL_THRESH), .N],
                merged_gene_tpm[(test_ID == tid) & (qval <= QVAL_THRESH) & (abs(log2fc) < LOG2FOLD_THRESH), .N],
                merged_gene_tpm[(test_ID == tid) & (qval <= QVAL_THRESH) & (abs(log2fc) >= LOG2FOLD_THRESH), .N]
            ),
            Total = rep(merged_gene_tpm[test_ID == tid, .N], 3)
        )
        dt[, Frac := N / Total]
        return(dt)
    }
))

counted_gene_tpm_plot <- dcast(
    counted_gene_tpm,
    test_ID ~ Significant + LargeFold,
    value.var = "Frac"
)

# ==============================================================================
# Plots
# ==============================================================================
gg_volcano_transcripts <- (
    ggplot(data = tested_transcripts)
    + geom_point(aes(x = b * log2(exp(1)), y = -log10(qval), colour = (qval < QVAL_THRESH & abs(b) >= log(2))))
    + geom_vline(aes(xintercept = -1), linetype = "dashed")
    + geom_vline(aes(xintercept = 1), linetype = "dashed")
    + geom_hline(aes(yintercept = -log10(QVAL_THRESH)), linetype = "dashed")
    + labs(x = expression(log[2] * "(Fold change)"), y = expression(-log[10] * "(FDR)"))
    + scale_colour_manual(
        breaks = c(FALSE, TRUE),
        values = c("#BDBDBD", "#000000")
    )
    + theme_minimal()
)
savefig(gg_volcano_transcripts, "Plots/volcano.transcripts")

gg_volcano_genes <- (
    ggplot(data = merged_gene_tpm)
    + geom_point(aes(
        x = log2fc,
        y = -log10(qval),
        colour = ((qval < QVAL_THRESH) & (abs(log2fc) >= 1))
    ))
    + geom_vline(aes(xintercept = -1), linetype = "dashed")
    + geom_vline(aes(xintercept = 1), linetype = "dashed")
    + geom_hline(aes(yintercept = -log10(QVAL_THRESH)), linetype = "dashed")
    + labs(x = expression(log[2] * "(Fold change)"), y = expression(-log[10] * "(FDR)"))
    + scale_colour_manual(
        breaks = c(FALSE, TRUE),
        values = c("#BDBDBD", "#000000")
    )
    + coord_cartesian(
        xlim = c(-7, 7),
        ylim = c(0, 20)
    )
    + theme_minimal()
)
savefig(gg_volcano_genes, "Plots/volcano.genes")

panel_height_ratios <- c(9, 1, 1)
gg_fc_bars <- ggplotGrob(
    ggplot(data = counted_gene_tpm)
    + geom_col(
        aes(
            x = factor(test_ID, levels = counted_gene_tpm_plot[order(FALSE_NA, TRUE_FALSE, TRUE_TRUE), test_ID], ordered = TRUE),
            y = Frac,
            fill = paste(Significant, LargeFold, sep = "_"),
        ),
        position = "stack"
    )
    + labs(x = NULL, y = "Genes in TAD around breakpoint (%)")
    + scale_fill_manual(
        breaks = c("FALSE_NA", "TRUE_FALSE", "TRUE_TRUE"),
        labels = c(
            "N.S.",
            expression("|" * log[2] * "(Fold Change)| < 1"),
            expression("|" * log[2] * "(Fold Change)| >= 1")
        ),
        values = c("#BDBDBD", "#DCA395", "#FF6347"),
        name = "Gene expression change"
    )
    + theme_minimal()
    + theme(
        axis.text.x = element_blank(),
        legend.position = "top",
        panel.grid.major.x = element_blank()
    )
)
gg_fc_ngenes <- ggplotGrob(
    ggplot(data = counted_gene_tpm[Significant == FALSE])
    + geom_col(aes(
        x = factor(test_ID, levels = counted_gene_tpm_plot[order(FALSE_NA, TRUE_FALSE, TRUE_TRUE), test_ID], ordered = TRUE),
        y = Total
    ))
    + labs(x = NULL)
    + scale_y_continuous(
        name = "Genes",
        limits = c(0, counted_gene_tpm[, max(Total)]),
        breaks = c(0, counted_gene_tpm[, max(Total)])
    )
    + theme_minimal()
    + theme(
        strip.text.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank()
    )
)
gg_fc_affecting_tad <- ggplotGrob(
    ggplot(data = tad_tests)
    + geom_col(aes(
        x = factor(test_ID, levels = counted_gene_tpm_plot[order(FALSE_NA, TRUE_FALSE, TRUE_TRUE), test_ID], ordered = TRUE),
        y = 1,
        fill = altered_TAD
    ))
    + labs(x = NULL)
    + scale_y_continuous(
        limits = c(0, 1),
        labels = NULL,
        name = "Altered\nTAD"
    )
    + scale_fill_manual(
        limits = c(TRUE, FALSE),
        labels = c("Yes", "No"),
        values = c("#000000", "#BDBDBD"),
        name = "SV affects TAD boundaries"
    )
    + guides(fill = FALSE)
    + theme_minimal()
    + theme(
        strip.text.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank()
    )
)
maxWidth <- grid::unit.pmax(
    gg_fc_bars$widths[2:5],
    gg_fc_ngenes$widths[2:5],
    gg_fc_affecting_tad$widths[2:5]
)
gg_fc_bars$widths[2:5] <- as.list(maxWidth)
gg_fc_ngenes$widths[2:5] <- as.list(maxWidth)
gg_fc_affecting_tad$widths[2:5] <- as.list(maxWidth)
gg_fc <- grid.arrange(
    gg_fc_bars,
    gg_fc_affecting_tad,
    gg_fc_ngenes,
    nrow = 3,
    heights = panel_height_ratios / sum(panel_height_ratios)
)
savefig(gg_fc, file.path(PLOT_DIR, "expression"))

# ==============================================================================
# Save tables
# ==============================================================================
fwrite(
    counted_gene_tpm,
    "1sample-results.tsv",
    sep = "\t",
    col.names = TRUE
)
