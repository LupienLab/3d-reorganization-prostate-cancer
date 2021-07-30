# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
source(file.path("..", "src", "savefig.R"))

PLOT_DIR <- "Plots"

# ==============================================================================
# Data
# ==============================================================================
# load sample metadata
metadata <- fread(
    "../../Data/External/LowC_Samples_Data_Available.tsv",
    sep = "\t",
    header = TRUE
)
metadata <- metadata[Include == "Yes"]
metadata[, SampleID := paste0("PCa", get("Sample ID"))]
SAMPLES <- metadata[, SampleID]


# load expression data for each tested gene
tested_genes <- fread(
    file.path("..", "2020-02-19_sv-disruption-expression", "sv-disruption-tests.expression.gene-level.tsv"),
    sep = "\t",
    header = TRUE
)


# ==============================================================================
# Analysis
# ==============================================================================
# ZNF516 case example
# --------------------------------------
znf516 <- melt(
    tested_genes[name == "ZNF516"],
    id.vars = "test_ID",
    measure.vars = SAMPLES,
    variable.name = "SampleID",
    value.name = "Expression"
)
znf516[, Mutated := FALSE]
znf516[SampleID == "PCa13848", Mutated := TRUE]

# add colour metadata
znf516 <- merge(
    x = znf516,
    y = metadata[, .SD, .SDcols = c("SampleID", "Patient ID", "Colour")],
    by = "SampleID"
)

# SLC4A4 case example
# --------------------------------------
slc4a4 <- melt(
    tested_genes[name == "SLC4A4"],
    id.vars = "test_ID",
    measure.vars = SAMPLES,
    variable.name = "SampleID",
    value.name = "Expression"
)
slc4a4[, Mutated := FALSE]
slc4a4[SampleID == "PCa53687", Mutated := TRUE]

# add colour metadata
slc4a4 <- merge(
    x = slc4a4,
    y = metadata[, .SD, .SDcols = c("SampleID", "Patient ID", "Colour")],
    by = "SampleID"
)

# ==============================================================================
# Plots
# ==============================================================================
# ZNF516
gg_znf516 <- (
    ggplot()
    + geom_boxplot(
        data = znf516[SampleID != "PCa13848"],
        aes(x = Mutated, y = log2(Expression)),
        alpha = 0.7,
        width = 0.3
    )
    + geom_point(
        data = znf516,
        aes(x = Mutated, y = log2(Expression), colour = SampleID),
        position = position_jitter(height = 0, width = 0.3)
    )
    + labs(y = expression(log[2] * " ZNF516 Expression (FPKM)"))
    + scale_x_discrete(
        breaks = c(FALSE, TRUE),
        labels = c("Others", "CPCG0366"),
        name = NULL
    )
    + scale_colour_manual(
        breaks = metadata[, SampleID],
        labels = metadata[, get("Patient ID")],
        values = metadata[, Colour],
    )
    + guides(colour = FALSE)
    + theme_minimal()
)
savefig(gg_znf516, file.path(PLOT_DIR, "ZNF516.expression"), width = 8)

# SLC4A4
gg_slc4a4 <- (
    ggplot()
    + geom_boxplot(
        data = slc4a4[SampleID != "PCa53687"],
        aes(x = Mutated, y = log2(Expression)),
        alpha = 0.7,
        width = 0.3
    )
    + geom_point(
        data = slc4a4,
        aes(x = Mutated, y = log2(Expression), colour = SampleID),
        position = position_jitter(height = 0, width = 0.3)
    )
    + labs(y = expression(log[2] * " SLC4A4 Expression (FPKM)"))
    + scale_x_discrete(
        breaks = c(FALSE, TRUE),
        labels = c("Others", "CPCG0339"),
        name = NULL
    )
    + scale_colour_manual(
        breaks = metadata[, SampleID],
        labels = metadata[, get("Patient ID")],
        values = metadata[, Colour],
    )
    + guides(colour = FALSE)
    + theme_minimal()
)
savefig(gg_slc4a4, file.path(PLOT_DIR, "SLC4A4.expression"), width = 8)
