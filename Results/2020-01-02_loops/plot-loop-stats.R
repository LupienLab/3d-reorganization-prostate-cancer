# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))

# ==============================================================================
# Data
# ==============================================================================
# load loop count table
loops = fread(file.path("Loops", "loop-stats.tsv"), sep = "\t", header = TRUE)

loops = loops[order(N_Loops), .SD]

# ==============================================================================
# Plots
# ==============================================================================
# set order for plotting
loops[, SampleID := factor(SampleID, levels = SampleID, ordered = TRUE)]
gg = (
    ggplot(data = loops)
    + geom_col(aes(x = SampleID, y = N_Loops))
    + labs(x = NULL, y = "Number of Loops")
    + theme_minimal()
    + theme(
        axis.text.x = element_text(angle = 90)
    )
)
ggsave(
    file.path("Plots", "loop-stats.png"),
    height = 12,
    width = 20,
    units = "cm"
)
