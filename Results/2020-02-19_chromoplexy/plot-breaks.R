# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("RCircos"))
suppressMessages(library("ggplot2"))
suppressMessages(library("scales"))

CHRS = paste0('chr', c(1:22, "X", "Y"))

hg38 = data.table(
    Chrom = factor(CHRS, levels = CHRS, ordered = TRUE),
    Length = c(
        248956422,
        242193529,
        198295559,
        190214555,
        181538259,
        170805979,
        159345973,
        145138636,
        138394717,
        133797422,
        135086622,
        133275309,
        114364328,
        107043718,
        101991189,
        90338345,
        83257441,
        80373285,
        58617616,
        64444167,
        46709983,
        50818468,
        156040895,
        57227415
    )
)
hg38[, Bins := ceiling(Length / 10^6)]

# ==============================================================================
# Data
# ==============================================================================
# load metadata
metadata <- fread("../../Data/External/LowC_Samples_Data_Available.tsv", sep = "\t", header = TRUE)
metadata[, SampleID := paste0("PCa", get("Sample ID"))]
SAMPLES <- metadata$SampleID

# load breakpoint data
breakpoints <- fread("Graphs/sv-breakpoints.tsv", sep = "\t", header = TRUE)
# convert chr to ordered factor
breakpoints[, chr := factor(chr, levels = CHRS, ordered = TRUE)]
breakpoint_pairs <- fread("Graphs/sv-breakpoints.paired.tsv", sep = "\t", header = TRUE)

# load hg38 ideogram
data(UCSC.HG38.Human.CytoBandIdeogram)

# ==============================================================================
# Analysis
# ==============================================================================
# convert breakpoints into Mbp bins and counts the rearrangements
breakpoints[, StartBin := floor(start / 10^6)]
breakpoints[, EndBin := floor(end / 10^6)]

# convert breakpoints in genomic regions into bins
breakpoints_binned = rbindlist(lapply(
    1:breakpoints[, .N],
    function(i) {
        start_bin = breakpoints[i, StartBin]
        end_bin = breakpoints[i, EndBin]
        dt = data.table(
            SampleID = breakpoints[i, SampleID],
            chr = breakpoints[i, chr],
            Bin = start_bin:end_bin,
            Count = 1
        )
        return(dt)
    }
))
# count all breakpoints in each bin
breakpoints_summed = breakpoints_binned[, sum(Count), by = c("SampleID", "chr", "Bin")]
colnames(breakpoints_summed) = c("SampleID", "chr", "Bin", "Count")

fwrite(
    breakpoints_summed[order(SampleID, chr, Bin)],
    "Statistics/breakpoints.binned.tsv",
    sep = "\t",
    col.names = TRUE
)

# count the number of inter-/intra-chromosomal pairs
inter_intra_counts <- breakpoint_pairs[, .N, by = "SampleID"]
inter_intra_counts$Intrachromosomal <- breakpoint_pairs[, chr_x == chr_y, by = "SampleID"][, sum(V1), by = "SampleID"]$V1
inter_intra_counts$Interchromosomal <- breakpoint_pairs[, chr_x != chr_y, by = "SampleID"][, sum(V1), by = "SampleID"]$V1

fwrite(
    inter_intra_counts,
    "Statistics/breakpoint-pairs.inter-intra-chromosomal.tsv",
    sep = "\t",
    col.names = TRUE
)

# calculate the length of events in each sample (i.e. size of each connected component)
breakpoint_components <- breakpoints[, .N, by = c("SampleID", "component_ID")]

# merge metadata for plotting
breakpoint_components <- merge(
    x = breakpoint_components,
    y = metadata[, .SD, .SDcols = c("SampleID", "T2E Status")],
    by = "SampleID"
)

breakpoint_components_counted <- breakpoint_components[, .N, by = c("SampleID", "T2E Status")]
colnames(breakpoint_components_counted)[3] <- "Events"
breakpoint_components_counted$Complex_Events <- breakpoint_components[,
    N > 2,
    by = "SampleID"
    ][, sum(V1), by = "SampleID"]$V1

# perform Mann-Whitney U test to see if the T2E samples have more complex events than the non-T2E samples
htest_events <- list(
    "total" = wilcox.test(
        x = breakpoint_components_counted[get("T2E Status") == "No", Events],
        y = breakpoint_components_counted[get("T2E Status") == "Yes", Events],
        alternative = "less"
    ),
    "complex" = wilcox.test(
        x = breakpoint_components_counted[get("T2E Status") == "No", Complex_Events],
        y = breakpoint_components_counted[get("T2E Status") == "Yes", Complex_Events],
        alternative = "less"
    )
)

fwrite(
    breakpoint_components_counted,
    "Statistics/breakpoint-components.tsv",
    sep = "\t",
    col.names = TRUE
)

# ==============================================================================
# Plots
# ==============================================================================
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
        paste0("Plots/circos/", s, ".circos.png"),
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
    # change colour of inter vs intra chromosomal events
    # params = RCircos.Get.Plot.Parameters()
    # params$PlotColor = breakpoints$Colour
    # RCircos.Reset.Plot.Parameters(params)
    dev.off()
    pdf(
        paste0("Plots/circos/", s, ".circos.pdf"),
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
    # change colour of inter vs intra chromosomal events
    # params = RCircos.Get.Plot.Parameters()
    # params$PlotColor = breakpoints$Colour
    # RCircos.Reset.Plot.Parameters(params)
    dev.off()
}

# 2. SV counts and summary statistics per patient
# --------------------------------------
# number of SVs detected per patient
gg_breakpoints = (
    ggplot(data = breakpoints[, .N, by = SampleID])
    + geom_col(aes(x = SampleID, y = N, fill = SampleID))
    + labs(x = NULL, y = "Unique Breakpoints")
    + guides(fill = FALSE)
    + theme_minimal()
    + theme(
        axis.text.x = element_text(angle = 90, vjust = 0, hjust = 0.5)
    )
)
ggsave(
    "Plots/breakpoint-stats/sv-counts.png",
    gg_breakpoints,
    height = 12,
    width = 20,
    units = "cm",
    dpi = 400
)
ggsave(
    "Plots/breakpoint-stats/sv-counts.pdf",
    gg_breakpoints,
    height = 12,
    width = 20,
    units = "cm",
    dpi = 400
)

gg_breakpoints_per_chrom = (
    ggplot(data = breakpoints[, .N, by = c("SampleID", "chr")])
    + geom_col(aes(x = SampleID, y = N, fill = SampleID))
    + labs(x = NULL, y = "Breakpoints")
    + guides(fill = guide_legend(title = "Patient"))
    + facet_grid(. ~ chr)
    + theme_minimal()
    + theme(
        panel.grid.major.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "bottom"
    )
)
ggsave(
    "Plots/breakpoint-stats/sv-counts-by-chrom.labelled.png",
    gg_breakpoints_per_chrom,
    height = 12,
    width = 40,
    units = "cm",
    dpi = 400
)
ggsave(
    "Plots/breakpoint-stats/sv-counts-by-chrom.labelled.pdf",
    gg_breakpoints_per_chrom,
    height = 12,
    width = 40,
    units = "cm",
    dpi = 400
)

# save as above without any text
gg_breakpoints_per_chrom_unlabelled = (
    gg_breakpoints_per_chrom 
    + labs(x = NULL, y = NULL)
    + guides(fill = guide_legend(title = NULL))
    + theme(
        axis.text.y = element_blank(),
        strip.text.x = element_blank()
    )
)
ggsave(
    "Plots/breakpoint-stats/sv-counts-by-chrom.png",
    gg_breakpoints_per_chrom_unlabelled,
    height = 12,
    width = 40,
    units = "cm",
    dpi = 400
)
ggsave(
    "Plots/breakpoint-stats/sv-counts-by-chrom.pdf",
    gg_breakpoints_per_chrom_unlabelled,
    height = 12,
    width = 40,
    units = "cm",
    dpi = 400
)

# number and location of SVs across all patients
gg_breakpoints_summed = (
    ggplot()
    + geom_col(
        data = breakpoints_summed,
        mapping = aes(x = Bin, y = Count, fill = SampleID),
        position = "stack"
    )
    # add geom_blank to ensure each facet is the correct size
    + geom_blank(data = hg38, mapping = aes(x = Bins, y = 0))
    + labs(x = "Position", y = "Count")
    + facet_grid(. ~ chr, scales = "free_x", space = "free_x", switch = "x")
    + theme_classic()
    + theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing = unit(0, "lines"),
        strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5)
    )
)
ggsave(
    "Plots/breakpoint-stats/sv-loci.png",
    gg_breakpoints_summed,
    height = 12,
    width = 40,
    units = "cm",
    dpi = 400
)
ggsave(
    "Plots/breakpoint-stats/sv-loci.pdf",
    gg_breakpoints_summed,
    height = 12,
    width = 40,
    units = "cm",
    dpi = 400
)

# megabase bins containing multiple breakpoints across samples
gg_recurrence = (
    ggplot(data = breakpoints_summed[, length(unique(SampleID)), by = c("chr", "Bin")])
    + geom_bar(aes(x = V1))
    + labs(x = "Number of patients", y = "Number of Mbp bins with a breakpoint")
    + scale_x_discrete(
        limits = seq(1, 13),
        breaks = seq(1, 13, by = 2),
        labels = seq(1, 13, by = 2)
    )
    + theme_minimal()
)
ggsave(
    "Plots/breakpoint-stats/sv-recurrence.png",
    gg_recurrence,
    height = 12,
    width = 20,
    units = "cm",
    dpi = 400
)
ggsave(
    "Plots/breakpoint-stats/sv-recurrence.pdf",
    gg_recurrence,
    height = 12,
    width = 20,
    units = "cm",
    dpi = 400
)

# 3. length and distribution of complex events
# --------------------------------------
gg_component_length = (
    ggplot(data = breakpoint_components)
    + geom_bar(
        aes(x = N, group = SampleID, fill = SampleID),
        colour = "#000000",
        position = "dodge"
    )
    + labs(x = "Breakpoints in SV", y = "Count")
    + guides(group = FALSE, fill = guide_legend(title = "Patients"))
    + theme_minimal()
)
ggsave(
    "plots/breakpoint-stats/sv-event-size.png",
    gg_component_length,
    height = 12,
    width = 20,
    units = "cm"
)
ggsave(
    "plots/breakpoint-stats/sv-event-size.pdf",
    gg_component_length,
    height = 12,
    width = 20,
    units = "cm"
)

# number of SV events detected in total
gg_t2e_components = (
    ggplot(data = breakpoint_components_counted)
    + geom_boxplot(aes(x = get("T2E Status"), y = Events, colour = get("T2E Status")), width = 0.5)
    + geom_point(
        aes(x = get("T2E Status"), y = Events, colour = get("T2E Status")),
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
        label = paste("p =", scientific(htest_events[["total"]]$p.value, digits = 3)),
        x = 1.5,
        y = 38,
        vjust = -1
    )
    + labs(x = NULL, y = "Structural variants")
    + guides(colour = FALSE)
    + scale_x_discrete(
        breaks = c("No", "Yes"),
        labels = c("ERG-", "ERG+")
    )
    + ylim(0, 40)
    + theme_minimal()
)
ggsave(
    "Plots/breakpoint-stats/sv-events.total.T2E-comparison.png",
    gg_t2e_components,
    height = 12,
    width = 20,
    units = "cm"
)
ggsave(
    "Plots/breakpoint-stats/sv-events.total.T2E-comparison.pdf",
    gg_t2e_components,
    height = 12,
    width = 20,
    units = "cm"
)

# number of complex components
gg_t2e_components_complex = (
    ggplot(data = breakpoint_components_counted)
    + geom_boxplot(aes(x = get("T2E Status"), y = Complex_Events, colour = get("T2E Status")), width = 0.5)
    + geom_point(
        aes(x = get("T2E Status"), y = Complex_Events, colour = get("T2E Status")),
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
        label = paste("p =", scientific(htest_events[["complex"]]$p.value, digits = 3)),
        x = 1.5,
        y = 13,
        vjust = -1
    )
    + labs(x = NULL, y = "Structural variants")
    + ylim(0, 15)
    + guides(colour = FALSE)
    + scale_x_discrete(
        breaks = c("No", "Yes"),
        labels = c("ERG-", "ERG+")
    )
    + theme_minimal()
)
ggsave(
    "Plots/breakpoint-stats/sv-events.complex.T2E-comparison.png",
    gg_t2e_components_complex,
    height = 12,
    width = 20,
    units = "cm"
)
ggsave(
    "Plots/breakpoint-stats/sv-events.complex.T2E-comparison.pdf",
    gg_t2e_components_complex,
    height = 12,
    width = 20,
    units = "cm"
)
