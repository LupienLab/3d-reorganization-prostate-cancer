# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("logging"))

loginfo("Loading packages")
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
split_comma_col <- function(v, f = identity) {
    # split into list
    splitv <- lapply(v, function(x) {
        strsplit(x, ",")[[1]]
    })
    # remove various non-informative characters (spaces, braces)
    splitv <- lapply(splitv, function(x) {
        gsub("[][ ]", "", x)
    })
    return(lapply(splitv, f))
}


# ==============================================================================
# Data
# ==============================================================================
loginfo("Loading data")

# load tests
sv_tests <- fread(
    "../2020-02-19_chromoplexy/Graphs/sv-disruption-tests.tsv",
    sep = "\t",
    header = TRUE
)
sv_tests$mut_samples <- split_comma_col(sv_tests$mut_samples)
sv_tests$nonmut_samples <- split_comma_col(sv_tests$nonmut_samples)

SAMPLES <- sv_tests[, sort(unique(unlist(mut_samples)))]

sv_bp_pairs <- fread(
    "../2020-02-19_chromoplexy/Graphs/sv-breakpoints.paired.tsv",
    sep = "\t",
    header = TRUE
)

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
            header = TRUE,
            fill = TRUE
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

transcript_table <- fread("all-samples.abundance.tsv", sep = "\t")


# ==============================================================================
# Analysis
# ==============================================================================
loginfo("Combining gene and transcript counts")

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
        gene_tpm <- transcript_table[target_id %in% transcript_ids, .(tpm = sum(tpm)), by = "sample"]
        gene_tpm[, gene_id := gene]

        # get test_ID
        tid <- as.integer(r["test_ID"])
        # add Mutated column for each patient involved in a test
        gene_tpm[, test_ID := tid]
        gene_tpm[, Mutated := ""]
        mut_samples <- sv_tests[test_ID == tid, unlist(mut_samples)]
        nonmut_samples <- sv_tests[test_ID == tid, unlist(nonmut_samples)]
        gene_tpm[sample %in% mut_samples, Mutated := "Mutated"]
        gene_tpm[sample %in% nonmut_samples, Mutated := "Nonmutated"]

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

# link tests together by the SV they come from
sv_events_tests <- sv_bp_pairs[,
    .(test_ID = c(test_ID_x, test_ID_y)),
    keyby = c("SampleID", "component_ID_x")
]

# PCa58215_3 and PCa58215_2 are actually part of the same event, so link them
# I'm changing this manually here as to not mess up the test_IDs from before.
# it was tough to get them in the first place and will cause a lot of headaches to reconstruct everything
sv_events_tests[(SampleID == "PCa58215") & (component_ID_x == 3), component_ID_x := 2]
sv_events_tests[, event_ID := paste(SampleID, component_ID_x, sep = "_")]

merged_gene_tpm_multi_IDs <- unique(merge(
    x = merged_gene_tpm,
    y = sv_events_tests[, .SD, .SDcols = c("event_ID", "test_ID")],
    by = "test_ID",
    allow.cartesian = TRUE
))

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

# classify the events as "only increased expression", "only decreased expression", "no changes", and "increased and descreased expression"
event_status <- rbindlist(lapply(
    unique(merged_gene_tpm_multi_IDs$event_ID),
    function(eid) {
        diff_expr_genes <- merged_gene_tpm_multi_IDs[event_ID == eid & qval < QVAL_THRESH]
        # fc_signs is then a vector with 0 (c()), 1 (c(1) or c(-1)), or 2 (c(-1, 1)) elements
        fc_signs <- sort(unique(diff_expr_genes[, sign(log2fc)]))
        status <- ifelse(
            length(fc_signs) == 0,
            "none",
            ifelse(
                length(fc_signs) == 2,
                "up and down",
                ifelse(fc_signs == -1, "down only", "up only")
            )
        )
        return(data.table(
            event_ID = eid,
            status = factor(status, ordered = TRUE, levels = c("up only", "up and down", "down only", "none"))
        ))
    }
))

merged_gene_tpm_multi_IDs <- merge(
    x = merged_gene_tpm_multi_IDs,
    y = event_status,
    by = "event_ID",
    all.x = TRUE
)

fwrite(merged_gene_tpm_multi_IDs, "summary-sv-disruption.tsv", sep = "\t")

# ==============================================================================
# Plots
# ==============================================================================
loginfo("Plotting figures")

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


# make two versions of this plot, one with labelled SVs and one without
gg_event_changes <- (
    ggplot(data = merged_gene_tpm_multi_IDs[status != "none"])
    + geom_point(
        aes(
            x = sign(log2fc) * pmin(-log10(qval), 16),
            y = event_ID,
            fill = qval < QVAL_THRESH
        ),
        shape = 21,
        size = 3
    )
    + geom_vline(
        xintercept = c(log10(QVAL_THRESH), -log10(QVAL_THRESH)),
        linetype = "dashed",
        colour = "#ff6347"
    )
    + labs(y = "Structural Variants")
    + coord_cartesian(
        xlim = c(-16, 16),
        expand = FALSE,
        clip = "off"
    )
    + scale_x_continuous(
        name = bquote(-log[10] * "(FDR) gene expression"),
        breaks = c(
            -16, -12, -8, -4,
            0,
            4, 8, 12, 16
        ),
        labels = c(
            "> 16", "12", "8", "4", "0", "4", "8", "12", "> 16"
        )
    )
    + scale_fill_manual(
        breaks = c(FALSE, TRUE),
        values = c("#bdbdbd", "#ff6347")
    )
    + facet_grid(status ~ ., scales = "free_y")
    + guides(fill = FALSE)
    + theme_minimal()
    + theme(
        panel.spacing = unit(2, "lines")
    )
)
# adjust facet heights
# adjusted from https://stackoverflow.com/questions/52341385/how-to-automatically-adjust-the-width-of-each-facet-for-facet-wrap
gp <- ggplotGrob(gg_event_changes)
facet_rows <- gp$layout$b[grepl("panel", gp$layout$name)]
# get unique y-axis vales per facet
y_events <- sapply(
    ggplot_build(gg_event_changes)$layout$panel_scales_y,
    function(l) length(l$range$range)
)
# change relative heights of facets
gp$heights[facet_rows] <- gp$heights[facet_rows] * y_events
savefig(gp, "Plots/event-changes.labelled")

# make the second version
gg_event_changes <- gg_event_changes + theme(axis.text.y = element_blank())
gp <- ggplotGrob(gg_event_changes)
facet_rows <- gp$layout$b[grepl("panel", gp$layout$name)]
y_events <- sapply(
    ggplot_build(gg_event_changes)$layout$panel_scales_y,
    function(l) length(l$range$range)
)
gp$heights[facet_rows] <- gp$heights[facet_rows] * y_events
savefig(gp, "Plots/event-changes")

# ==============================================================================
# Save tables
# ==============================================================================
fwrite(
    counted_gene_tpm,
    "1sample-results.tsv",
    sep = "\t",
    col.names = TRUE
)
