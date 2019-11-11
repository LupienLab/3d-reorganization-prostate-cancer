# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
suppressMessages(library("grid"))

CHRS = paste0("chr", c(1:22, "X", "Y"))
CHRS = factor(CHRS, levels = CHRS, ordered = TRUE)

# ==============================================================================
# Data
# ==============================================================================
# read minimum pairwise symmetric similarity
mpss = fread("TAD-comparisons/w_30.40000bp.similarity.tsv")
mpss[, chr_nat := factor(chr, levels = CHRS)]

# ==============================================================================
# Plots
# ==============================================================================
gg = (
    ggplot(data = mpss)
    + geom_rect(aes(xmin = start, ymin = 0, xmax = end, ymax = similarity))
    + labs(x = "Position", y = "Min. Sym. Pairwise Similarity")
    + facet_grid(. ~ chr_nat, switch = "x", scales = "free_x", space = "free")
    + theme_minimal()
    + theme(
        panel.spacing = unit(0, "lines"),
        panel.border = element_rect(size = 0.1, colour = "black", fill = NA),
        strip.text.x = element_text(angle=90, hjust=0.5, vjust=0.5),
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
    )
)

ggsave(
    "Plots/tad-similarity-track.png",
    height = 12,
    width = 20,
    units = "cm"
)
