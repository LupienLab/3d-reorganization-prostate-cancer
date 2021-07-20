# ==============================================================================
# Meta
# ==============================================================================
# Plot gene expression
# --------------------------------------
# Description: Plot the gene expression changes for the T2E translocation insertion site
# Author: James Hawley


# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("logging"))

loginfo("Loading packages")
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
suppressMessages(library("sleuth"))
source("helper-functions.R")

RES_DIR <- file.path(
    "..", "..", "results", "2020-06-18_sv-disruption-expression"
)
PLOT_DIR <- file.path(RES_DIR, "Plots")


# ==============================================================================
# Data
# ==============================================================================
loginfo("Loading data")

# load Kallisto/Sleuth object
so <- list(
    "ERG" = file.path(RES_DIR, "sleuth", "test_46.sleuth-object.rds"),
    "TMPRSS2" = file.path(RES_DIR, "sleuth", "test_45.sleuth-object.rds"),
    "ZNF516" = file.path(RES_DIR, "sleuth", "test_42.sleuth-object.rds"),
    "PMEPA1" = file.path(RES_DIR, "sleuth", "test_43.sleuth-object.rds")
)

mut_samples <- list(
    "ERG" = c("PCa13848", "PCa19121", "PCa3023", "PCa40507", "PCa51852", "PCa56413"),
    "TMPRSS2" = c("PCa13848", "PCa19121", "PCa3023", "PCa40507", "PCa51852", "PCa56413"),
    "ZNF516" = c("PCa13848"),
    "PMEPA1" = c("PCa13848")
)
nonmut_samples <- lapply(
    names(mut_samples),
    function(alt_gene) {
        return(setdiff(
            c(
                "PCa13266", "PCa13848", "PCa14121", "PCa19121", "PCa3023", "PCa33173",
                "PCa40507", "PCa51852", "PCa53687", "PCa56413", "PCa57294", "PCa58215"
            ),
            mut_samples[[alt_gene]]
        ))
    }
)
names(nonmut_samples) <- names(mut_samples)

# load GENCODE annotations
gencode_genes <- fread(
    file.path("..", "..", "data", "External", "GENCODE", "gencode.v33.all-genes.bed"),
    sep = "\t",
    header = FALSE,
    col.names = c("chr", "start", "end", "strand", "ens_gene", "gene_name")
)
gencode_tx <- fread(
    file.path("..", "..", "data", "External", "GENCODE", "gencode.v33.all-transcripts.bed"),
    sep = "\t",
    header = FALSE,
    col.names = c("chr", "start", "end", "strand", "ens_gene", "gene_name", "ens_transcript", "transcript_name")
)

# ==============================================================================
# Analysis
# ==============================================================================
loginfo("Performing calculations")

ALT_GENES <- names(so)

# calculate aggregated gene count
so_counts <- lapply(
    ALT_GENES,
    function(alt_gene) {
        get_important_sleuth_info(
            path = so[[alt_gene]],
            alt_gene = alt_gene,
            base = 2
        )
    }
)
names(so_counts) <- ALT_GENES

# ==============================================================================
# Plots
# ==============================================================================
loginfo("Plotting data")

gene_pairs <- list(
    c("ERG", "TMPRSS2"),
    c("ZNF516", "PMEPA1")
)

for (gp in gene_pairs) {
    # get the estimated gene counts for each gene in the pair
    quant <- unique(rbindlist(lapply(gp, function(alt_gene) so_counts[[alt_gene]]$gene_counts)))
    dge <- unique(rbindlist(lapply(gp, function(alt_gene) so_counts[[alt_gene]]$gene_de)))
    # save this particular table for posterity
    fwrite(
        quant,
        file.path(
            RES_DIR,
            paste0("quant.", paste(gp, collapse= "-"), ".tsv")
        ),
        sep = "\t",
        col.names = TRUE
    )
    fwrite(
        dge,
        file.path(
            RES_DIR,
            paste0("test.", paste(gp, collapse= "-"), ".tsv")
        ),
        sep = "\t",
        col.names = TRUE
    )
    # combined mut samples
    ms <- unique(unlist(mut_samples[gp]))
    nms <- unique(unlist(nonmut_samples[gp]))
    gg <- (
        ggplot()
        + geom_boxplot(
            data = quant[sample %in% nms],
            mapping = aes(x = gene_name, y = est_counts),
            outlier.shape = NA
        )
        + geom_point(
            data = quant[sample %in% nms],
            mapping = aes(x = gene_name, y = est_counts),
            fill = "#BDBDBD",
            colour = "#000000",
            shape = 21,
            size = 4,
            position = position_jitter(height = 0, width = 0.3)
        )
        + geom_point(
            data = quant[sample %in% ms],
            mapping = aes(x = gene_name, y = est_counts),
            fill = "#ff6347",
            colour = "#000000",
            shape = 21,
            size = 4,
            position = position_jitter(height = 0, width = 0.3)
        )
        + labs(
            x = "Gene",
            y = expression("Estimated mRNA count (" * log[2] * ")")
        )
        + theme_minimal()
        + jrh_theme()
    )
    savefig(
        gg,
        file.path(
            PLOT_DIR, 
            paste0("gene-exprs.", paste(gp, collapse = "-"))
        ),
        height = 10,
        width = 6
    )
}
