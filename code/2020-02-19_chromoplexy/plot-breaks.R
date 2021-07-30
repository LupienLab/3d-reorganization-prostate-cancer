# ==============================================================================
# Meta
# ==============================================================================
# plot-breaks
# ------------------------------------------------
# Author: James Hawley
# Description: Plot SV breakpoints identified in Hi-C data


# ==============================================================================
# Environment
# ==============================================================================
suppressWarnings(library("logging"))
loginfo("Loading packages")

suppressWarnings(library("data.table"))
suppressWarnings(library("ggplot2"))
suppressMessages(library("RCircos"))
suppressMessages(library("scales"))
source(file.path("..", "src", "savefig.R"))

RES_DIR <- file.path("..", "..", "results", "2020-02-19_chromoplexy")
GRAPH_DIR <- file.path(RES_DIR, "Graphs")
STATS_DIR <- file.path(RES_DIR, "Statistics")
PLOT_DIR <- file.path(RES_DIR, "Plots")

# ==============================================================================
# Data
# ==============================================================================
loginfo("Loading data")

CHRS <- paste0("chr", c(1:22, "X", "Y"))
hg38 <- fread(
    file.path(
        "..", "..", "data", "Processed", "2019-06-18_PCa-LowC-sequencing",
        "hg38.sizes.txt"
    ),
    sep = "\t",
    header = FALSE,
    col.names = c("Chrom", "Length")
)
hg38[, Bins := ceiling(Length / 10^6)]

# load metadata
metadata <- fread(
    file.path(
        "..", "..", "data", "External", "LowC_Samples_Data_Available.tsv"
    ),
    sep = "\t",
    header = TRUE
)
metadata <- metadata[Include == "Yes"]
metadata[, SampleID := paste0("PCa", get("Sample ID"))]
SAMPLES <- metadata$SampleID

# load breakpoint data
breakpoints <- fread(
    file.path(GRAPH_DIR, "sv-breakpoints.tsv"),
    sep = "\t",
    header = TRUE
)
# convert chr to ordered factor
breakpoints[, chr := factor(chr, levels = CHRS, ordered = TRUE)]

breakpoint_pairs <- fread(
    file.path(GRAPH_DIR, "sv-breakpoints.paired.tsv"),
    sep = "\t",
    header = TRUE
)

# load hg38 ideogram
data(UCSC.HG38.Human.CytoBandIdeogram)

# load summary data
breakpoints_by_chrom <- fread(
    file.path(STATS_DIR, "breakpoints.by-chrom.tsv"),
    sep = "\t",
    header = TRUE,
    col.names = c("SampleID", "T2E_Status", "chr", "N", "N_per_mb")
)

breakpoints_summed <- fread(
    file.path(STATS_DIR, "breakpoints.binned.tsv")
)

inter_intra_counts <- fread(
    file.path(STATS_DIR, "breakpoint-pairs.inter-intra-chromosomal.tsv")
)

sv_components <- fread(
    file.path(STATS_DIR, "sv-components.tsv")
)

sv_components_counted <- fread(
    file.path(STATS_DIR, "sv-components.counts.tsv")
)

htests <- readRDS(file.path(STATS_DIR, "htests.rds"))

# ==============================================================================
# Analysis
# ==============================================================================
loginfo("Calculating breakpoint pair distances")

# calculate the distance between intra-chromosomal breakpoint pairs
breakpoint_pairs[
    chr_x == chr_y,
    dist := pmin(
        abs(start_y - end_x),
        abs(start_x - end_y)
    )
]

# ==============================================================================
# Plots
# ==============================================================================
loginfo("Plotting figures")
# 1. Circos plots
# --------------------------------------
# set core component for circos plots
RCircos.Set.Core.Components(
    UCSC.HG38.Human.CytoBandIdeogram,
    NULL,
    tracks.inside = 0,
    tracks.outside = 0
)
# plot circos plots for each patient
for (s in SAMPLES) {
    cat(s, "\n")
    png(
        file.path(PLOT_DIR, "circos", paste0(s, ".circos.png")),
        width = 12,
        height = 12,
        units = "cm",
        res = 400
    )
    RCircos.Set.Plot.Area()
    RCircos.Chromosome.Ideogram.Plot()
    RCircos.Link.Plot(
        link.data = breakpoint_pairs[
            SampleID == s,
            .(chr_x, start_x, end_x, chr_y, start_y, end_y)
        ],
        track.num = 1,
        by.chromosome = TRUE
    )
    dev.off()
    pdf(
        file.path(PLOT_DIR, "circos", paste0(s, ".circos.pdf")),
        width = 12,
        height = 12
    )
    RCircos.Set.Plot.Area()
    RCircos.Chromosome.Ideogram.Plot()
    RCircos.Link.Plot(
        link.data = breakpoint_pairs[
            SampleID == s,
            .(chr_x, start_x, end_x, chr_y, start_y, end_y)
        ],
        track.num = 1,
        by.chromosome = TRUE
    )
    dev.off()
}

# 2. SV counts and summary statistics per patient
# --------------------------------------
# number of SVs detected per patient
gg_breakpoints <- (
    ggplot(data = breakpoints[, .N, by = SampleID])
    + geom_col(aes(x = SampleID, y = N, fill = SampleID))
    + labs(x = NULL, y = "Unique Breakpoints")
    + scale_x_discrete(
        breaks = metadata[, SampleID],
        labels = metadata[, get("Patient ID")]
    )
    + scale_y_continuous(
        limits = c(0, 100)
    )
    + scale_fill_manual(
        breaks = metadata[, SampleID],
        labels = metadata[, get("Patient ID")],
        values = metadata[, Colour],
        name = "Patient"
    )
    + guides(fill = FALSE)
    + theme_minimal()
    + theme(
        axis.text.x = element_text(angle = 90, vjust = 0, hjust = 0.5)
    )
)
savefig(
    gg_breakpoints,
    file.path(PLOT_DIR, "breakpoint-stats", "breakpoint-counts")
)

gg_breakpoints_per_chrom <- (
    ggplot(data = breakpoints_by_chrom)
    + geom_col(
        aes(x = chr, y = N_per_mb, fill = T2E_Status),
        position = "stack",
        colour = "#000000"
    )
    + labs(x = NULL, y = "Breakpoints / Mb")
    + scale_x_discrete(
        breaks = CHRS,
        labels = CHRS,
        drop = FALSE
    )
    + scale_fill_manual(
        breaks = c("No", "Yes"),
        labels = c("T2E-", "T2E+"),
        values = c("#418B3D", "#3215C1"),
        name = "T2E Status"
    )
    + theme_minimal()
    + theme(
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "bottom"
    )
)
savefig(
    gg_breakpoints_per_chrom,
    file.path(
        PLOT_DIR,
        "breakpoint-stats",
        "breakpoint-counts.by-chrom.labelled"
    )
)

# save as above without any text
gg_breakpoints_per_chrom_unlabelled <- (
    gg_breakpoints_per_chrom
    + labs(x = NULL, y = NULL)
        + guides(fill = guide_legend(title = NULL))
        + theme(
            axis.text.y = element_blank(),
            strip.text.x = element_blank()
        )
)
savefig(
    gg_breakpoints_per_chrom_unlabelled,
    file.path(
        PLOT_DIR,
        "breakpoint-stats",
        "breakpoint-counts.by-chrom"
    )
)

# distance distribution between intra-chromosomal breakpoint pairs
gg_breakpoint_pair_dist <- (
    ggplot(data = breakpoint_pairs)
    + geom_histogram(aes(x = dist / 1e6, binwidth = 10))
    + geom_density(aes(x = dist / 1e6, y = 10 * ..scaled..))
    + labs(x = "Breakpoint pair distance (Mbp)", y = "Frequency")
    + theme_minimal()
)
savefig(
    gg_breakpoint_pair_dist,
    file.path(
        PLOT_DIR,
        "breakpoint-stats",
        "breakpoint-pairs.distance-distribution"
    )
)


# megabase bins containing multiple breakpoints across samples
gg_recurrence <- (
    ggplot(
        data = breakpoints_summed[,
            length(unique(SampleID)),
            by = c("chr", "Bin")
        ]
    )
    + geom_bar(aes(x = V1))
    + labs(
        x = "# Patients sharing a mutated Mbp bin",
        y = "# Mbp bins with a recurrent breakpoint"
    )
    + scale_x_discrete(
        limits = seq(1, 6),
        breaks = seq(1, 6),
        labels = seq(1, 6)
    )
    + theme_minimal()
)
savefig(
    gg_recurrence,
    file.path(
        PLOT_DIR,
        "breakpoint-stats",
        "breakpoint-recurrence"
    )
)

# total number of inter-/intra-chromosomal breakpoint pairs
gg_inter_intra <- (
    ggplot(data = inter_intra_counts[Class != "Total"])
    + geom_col(
        aes(x = SampleID, y = Count, fill = Class),
        colour = "#000000",
        position = "dodge"
    )
    + scale_fill_manual(
        breaks = c("Interchromosomal", "Intrachromosomal"),
        labels = c("Inter-chromosomal", "Intra-chromosomal"),
        values = c("#3F3FFF", "#FF7F7F"),
        name = ""
    )
    + scale_x_discrete(
        breaks = metadata[, SampleID],
        labels = metadata[, get("Patient ID")]
    )
    + labs(x = NULL, y = "Breakpoints")
    + theme_minimal()
    + theme(
        axis.text.x = element_text(angle = 90),
        legend.position = "bottom"
    )
)
savefig(
    gg_inter_intra,
    file.path(
        PLOT_DIR,
        "breakpoint-stats",
        "breakpoint-pairs.inter-intra-chromosomal"
    )
)

gg_inter_intra_t2e <- (
    ggplot(
        data = inter_intra_counts[
            Class != "Total",
            .(N = sum(Count)),
            by = c("T2E_Status", "Class")
        ]
    )
    + geom_col(
        aes(x = T2E_Status, y = N, fill = Class),
        colour = "#000000",
        position = "dodge"
    )
    + scale_x_discrete(
        breaks = c("Yes", "No"),
        labels = c("T2E+", "T2E-")
    )
    + scale_fill_manual(
        breaks = c("Interchromosomal", "Intrachromosomal"),
        labels = c("Inter-chromosomal", "Intra-chromosomal"),
        values = c("#3F3FFF", "#FF7F7F"),
        name = ""
    )
    + labs(x = NULL, y = "Total Breakpoint Pairs")
    + coord_flip()
    + theme_minimal()
    + theme(
        axis.text.x = element_text(angle = 90),
        legend.position = "bottom"
    )
)
savefig(
    gg_inter_intra_t2e,
    file.path(
        PLOT_DIR,
        "breakpoint-stats",
        "breakpoint-pairs.inter-intra-chromosomal.T2E-total"
    )
)

gg_inter_intra_t2e_comp <- (
    ggplot(data = inter_intra_counts[Class != "Total"])
    + geom_boxplot(
        aes(x = T2E_Status, y = Count, colour = T2E_Status),
        alpha = 0.2,
        outlier.shape = NA
    )
    + geom_point(
        aes(x = T2E_Status, y = Count, colour = T2E_Status),
        position = position_jitter(height = 0, width = 0.4)
    )
    + scale_x_discrete(
        breaks = c("Yes", "No"),
        labels = c("T2E+", "T2E-")
    )
    + scale_colour_manual(
        breaks = c("No", "Yes"),
        labels = c("T2E-", "T2E+"),
        values = c("#418B3D", "#3215C1"),
        name = ""
    )
    + labs(x = NULL, y = "Breakpoint Pairs")
    + guides(colour = FALSE)
    + coord_flip()
    + facet_wrap(~Class, ncol = 1)
    + theme_minimal()
    + theme(
        axis.text.x = element_text(angle = 90),
        legend.position = "bottom"
    )
)
savefig(
    gg_inter_intra_t2e_comp,
    file.path(
        PLOT_DIR,
        "breakpoint-stats",
        "breakpoint-pairs.inter-intra-chromosomal.T2E-comparison"
    )
)


# 3. length and distribution of complex events
# --------------------------------------
gg_component_length <- (
    ggplot(data = sv_components)
    + geom_bar(
        aes(x = N_Breakpoints, group = SampleID, fill = SampleID),
        position = position_dodge(preserve = "single")
    )
    + scale_fill_manual(
        breaks = metadata[, SampleID],
        labels = metadata[, get("Patient ID")],
        values = metadata[, Colour],
        name = "Patient"
    )
    + labs(x = "Breakpoints in SV", y = "Count")
    + guides(group = FALSE, fill = guide_legend(title = "Patient"))
    + theme_minimal()
)
savefig(
    gg_component_length,
    file.path(
        PLOT_DIR,
        "breakpoint-stats",
        "sv-events.size"
    )
)

# number of SV events detected in total
gg_t2e_components <- (
    ggplot(data = sv_components_counted)
    + geom_boxplot(
        aes(x = get("T2E Status"), y = N_Events, colour = get("T2E Status")),
        width = 0.5,
        outlier.shape = NA
    )
    + geom_point(
        aes(x = get("T2E Status"), y = N_Events, colour = get("T2E Status")),
        position = position_jitter(width = 0.1, height = 0)
    )
    + geom_path(
        data = data.table(
            x = rep(c("No", "Yes"), each = 2),
            y = c(10, 38, 38, 37)
        ),
        aes(x = x, y = y, group = 1)
    )
    + annotate(
        geom = "text",
        label = paste("p =", scientific(htests[["total"]]$p.value, digits = 3)),
        x = 1.5,
        y = 38,
        vjust = -1
    )
    + labs(x = NULL, y = "Structural variants")
    + scale_x_discrete(
        breaks = c("No", "Yes"),
        labels = c("T2E-", "T2E+")
    )
    + scale_colour_manual(
        breaks = c("No", "Yes"),
        labels = c("T2E-", "T2E+"),
        values = c("#418B3D", "#3215C1")
    )
    + guides(colour = FALSE)
    + ylim(0, 40)
    + theme_minimal()
)
savefig(
    gg_t2e_components,
    file.path(
        PLOT_DIR,
        "breakpoint-stats",
        "sv-events.total.T2E-comparison"
    ),
    height = 20,
    width = 12
)

# number of complex components
gg_t2e_components_complex <- (
    ggplot(data = sv_components_counted)
    + geom_boxplot(
        aes(
            x = get("T2E Status"),
            y = N_Complex_Events,
            colour = get("T2E Status")
        ),
        width = 0.5
    )
    + geom_point(
        aes(
            x = get("T2E Status"),
            y = N_Complex_Events,
            colour = get("T2E Status")
        ),
        position = position_jitter(width = 0.1, height = 0)
    )
    + geom_path(
        data = data.table(
            x = rep(c("No", "Yes"), each = 2),
            y = c(3, 13, 13, 12)
        ),
        aes(x = x, y = y, group = 1)
    )
    + annotate(
        geom = "text",
        label = paste("p =", scientific(htests[["complex"]]$p.value, digits = 3)),
        x = 1.5,
        y = 13,
        vjust = -1
    )
    + labs(x = NULL, y = "Structural variants")
    + ylim(0, 15)
    + guides(colour = FALSE)
    + scale_x_discrete(
        breaks = c("No", "Yes"),
        labels = c("T2E-", "T2E+")
    )
    + scale_colour_manual(
        breaks = c("No", "Yes"),
        labels = c("T2E-", "T2E+"),
        values = c("#2A363B", "#019875")
    )
    + theme_minimal()
)
savefig(
    gg_t2e_components_complex,
    file.path(
        PLOT_DIR,
        "breakpoint-stats",
        "sv-events.complex.T2E-comparison"
    ),
    height = 20,
    width = 12
)

# number of chromosomes that events span
gg_sv_chrom_span <- (
    ggplot(data = sv_components[, .N, by = N_Chr])
    + geom_col(aes(x = N_Chr, y = N, fill = N_Chr))
    + labs(x = "Chromosomes in a complex SV", y = "Count")
    + guides(fill = FALSE)
    + scale_fill_viridis_c()
    + theme_minimal()
)
savefig(
    gg_sv_chrom_span,
    file.path(
        PLOT_DIR,
        "breakpoint-stats",
        "sv-events.chromosomes"
    )
)
