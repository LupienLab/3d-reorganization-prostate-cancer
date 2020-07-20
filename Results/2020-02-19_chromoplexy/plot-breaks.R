# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("RCircos"))
suppressMessages(library("ggplot2"))
suppressMessages(library("scales"))
source("plotting-helper.R")

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
metadata <- metadata[Include == "Yes"]
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

# add T2E status
breakpoints <- merge(
    breakpoints,
    metadata[, .SD, .SDcols = c("SampleID", "T2E Status")],
    by = "SampleID"
)

# count breakpoints per sample
breakpoints_counted <- breakpoints[, .N, keyby = c("SampleID", "T2E Status")]

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

breakpoints_by_chrom = breakpoints[, .N, keyby = c("SampleID", "chr")]
breakpoints_by_chrom[, N_per_mb := apply(.SD, 1, function(r) {as.numeric(r["N"]) / hg38[Chrom == r["chr"], (Length / 10^6)]})]

# merge T2E information
breakpoints_by_chrom <- merge(
    x = breakpoints_by_chrom,
    y = metadata[, .(SampleID, T2E_Status = get("T2E Status"))],
    by = "SampleID"
)
fwrite(
    breakpoints_by_chrom[order(SampleID, chr)],
    "Statistics/breakpoints.by-chrom.tsv",
    sep = "\t",
    col.names = TRUE
)


# count the number of inter-/intra-chromosomal pairs
inter_intra_counts <- breakpoint_pairs[, .(Total = .N), by = "SampleID"]
inter_intra_counts$Intrachromosomal <- breakpoint_pairs[, chr_x == chr_y, by = "SampleID"][, sum(V1), by = "SampleID"]$V1
inter_intra_counts$Interchromosomal <- breakpoint_pairs[, chr_x != chr_y, by = "SampleID"][, sum(V1), by = "SampleID"]$V1
inter_intra_counts$T2E_Status <- metadata[, get("T2E Status")] # these are in the same order, so I don't need to merge
inter_intra_counts <- melt(
    inter_intra_counts,
    id.vars = c("SampleID", "T2E_Status"),
    variable.name = "Class",
    value.name = "Count"
)

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
htest_breaks <- wilcox.test(
    x = breakpoints_counted[get("T2E Status") == "No", N],
    y = breakpoints_counted[get("T2E Status") == "Yes", N],
    alternative = "less"
)

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

# calculate the number of chromosomes involved in an SV
sv_chrom_span <- breakpoints[, length(unique(chr)), by = c("SampleID", "component_ID")]
sv_chrom_span <- merge(
    sv_chrom_span,
    breakpoint_components,
    by = c("SampleID", "component_ID")
)


fisher.test(
    rbind(
        c(inter_intra_counts[Class == "Intrachromosomal" & T2E_Status == "No", sum(Count)], inter_intra_counts[Class == "Intrachromosomal" & T2E_Status == "Yes", sum(Count)]),
        c(inter_intra_counts[Class == "Interchromosomal" & T2E_Status == "No", sum(Count)], inter_intra_counts[Class == "Interchromosomal" & T2E_Status == "Yes", sum(Count)])
    ),
    alternative = "two.sided"
)

breakpoint_pairs[chr_x == chr_y, dist := pmin(abs(start_y - end_x), abs(start_x - end_y))]

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
savefig(gg_breakpoints, "Plots/breakpoint-stats/sv-counts")

gg_breakpoints_per_chrom = (
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
savefig(gg_breakpoints_per_chrom, "Plots/breakpoint-stats/sv-counts.by-chrom.labelled")

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
savefig(gg_breakpoints_per_chrom_unlabelled, "Plots/breakpoint-stats/sv-counts.by-chrom")

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
    + scale_fill_manual(
        breaks = metadata[, SampleID],
        labels = metadata[, get("Patient ID")],
        values = metadata[, Colour],
        name = "Patient"
    )
    + facet_grid(. ~ chr, scales = "free_x", space = "free_x", switch = "x")
    + theme_classic()
    + theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing = unit(0, "lines"),
        strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5)
    )
)
savefig(gg_breakpoints_summed, "Plots/breakpoint-stats/sv-loci", width = 40)

gg_breakpoint_pair_dist <- (
    ggplot(data = breakpoint_pairs)
    + geom_histogram(aes(x = dist / 1e6, binwidth = 10))
    + geom_density(aes(x = dist / 1e6, y = 10 * ..scaled..))
    + labs(x = "Breakpoint pair distance (Mbp)", y = "Frequency")
    + theme_minimal()
)
savefig(gg_breakpoint_pair_dist, "Plots/breakpoint-stats/sv-pair-distance-distribution")


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
savefig(gg_recurrence, "Plots/breakpoint-stats/sv-recurrence")

gg_inter_intra = (
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
savefig(gg_inter_intra, "Plots/breakpoint-stats/breakpoint-pairs.inter-intra-chromosomal")

gg_inter_intra_t2e = (
    ggplot(data = inter_intra_counts[Class != "Total", .(N = sum(Count)), by = c("T2E_Status", "Class")])
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
savefig(gg_inter_intra_t2e, "Plots/breakpoint-stats/breakpoint-pairs.inter-intra-chromosomal.T2E-total")

gg_inter_intra_t2e_comp = (
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
    + facet_wrap(~ Class, ncol = 1)
    + theme_minimal()
    + theme(
        axis.text.x = element_text(angle = 90),
        legend.position = "bottom"
    )
)
savefig(gg_inter_intra_t2e_comp, "Plots/breakpoint-stats/breakpoint-pairs.inter-intra-chromosomal.T2E-comparison")


# 3. length and distribution of complex events
# --------------------------------------
gg_component_length = (
    ggplot(data = breakpoint_components)
    + geom_bar(
        aes(x = N, group = SampleID, fill = SampleID),
        position = position_dodge(preserve = "single")
    )
    + scale_fill_manual(
        breaks = metadata[, SampleID],
        labels = metadata[, get("Patient ID")],
        values = metadata[, Colour],
        name = "Patient"
    )
    + labs(x = "Breakpoints in SV", y = "Count")
    + guides(group = FALSE, fill = guide_legend(title = "Patients"))
    + theme_minimal()
)
savefig(gg_component_length, "Plots/breakpoint-stats/sv-events.size")

# number of SV events detected in total
gg_t2e_components = (
    ggplot(data = breakpoint_components_counted)
    + geom_boxplot(
        aes(x = get("T2E Status"), y = Events, colour = get("T2E Status")),
        width = 0.5,
        outlier.shape = NA
    )
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
savefig(gg_t2e_components, "Plots/breakpoint-stats/sv-events.total.T2E-comparison", height = 20, width = 12)

gg_t2e_components_bps = (
    ggplot(data = breakpoints[, .N, by = c("SampleID", "T2E Status")])
    + geom_boxplot(
        aes(x = get("T2E Status"), y = N, colour = get("T2E Status")),
        width = 0.5,
        outlier.shape = NA
    )
    + geom_point(
        aes(x = get("T2E Status"), y = N, colour = get("T2E Status")),
        position = position_jitter(width = 0.1, height = 0)
    )
    + geom_path(
        data = data.table(
            x = rep(c("No", "Yes"), each = 2),
            y = c(25, 98, 98, 96)
        ),
        aes(x = x, y = y, group = 1)
    )
    + annotate(
        geom = "text",
        label = paste("p =", scientific(htest_breaks$p.value, digits = 3)),
        x = 1.5,
        y = 98,
        vjust = -1
    )
    + labs(x = NULL, y = "Unique Breakpoints")
    + scale_x_discrete(
        breaks = c("No", "Yes"),
        labels = c("T2E-", "T2E+")
    )
    + scale_y_continuous(
        limits = c(0, 100)
    )
    + scale_colour_manual(
        breaks = c("No", "Yes"),
        labels = c("T2E-", "T2E+"),
        values = c("#2A363B", "#019875")
    )
    + guides(colour = FALSE)
    + theme_minimal()
)
savefig(gg_t2e_components_bps, "Plots/breakpoint-stats/sv-breakpoints.T2E-comparison", height = 20, width = 12)

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
        labels = c("T2E-", "T2E+")
    )
    + scale_colour_manual(
        breaks = c("No", "Yes"),
        labels = c("T2E-", "T2E+"),
        values = c("#2A363B", "#019875")
    )
    + theme_minimal()
)
savefig(gg_t2e_components_complex, "Plots/breakpoint-stats/sv-events.complex.T2E-comparison", height = 20, width = 12)

# number of chromosomes that events span
gg_sv_chrom_span <- (
    ggplot(data = sv_chrom_span[N > 2, .(V2 = .N), keyby = "V1"])
    + geom_col(aes(x = V1, y = V2, fill = V1))
    + labs(x = "Chromosomes in a complex SV", y = "Count")
    + guides(fill = FALSE)
    + scale_fill_viridis_c()
    + theme_minimal()
)
savefig(gg_sv_chrom_span, "Plots/breakpoint-stats/sv-events.chromosomes")
