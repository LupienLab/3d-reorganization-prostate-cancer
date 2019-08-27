# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
suppressMessages(library("argparse"))
suppressMessages(library("RCircos"))
suppressMessages(library("ggbio"))

if (!interactive()) {
    PARSER <- argparse::ArgumentParser(
        description = "Make circos and other plots for breakpoint calls"
    )
    PARSER$add_argument(
        "breaks",
        type = "character",
        help = "Aggregated breakpoint call TSV file"
    )
    PARSER$add_argument(
        "-o", "--output",
        type = "character",
        help = "Output summed breakpoints file",
        default = "Breakpoints/Default/breakpoints.binned.tsv"
    )
    ARGS <- PARSER$parse_args()
} else {
    ARGS = list(
        breaks = "Breakpoints/Default/breakpoints.tsv",
        output = "Breakpoints/Default/breakpoints.binned.tsv"
    )
}

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
hg38[, Bins := round(Length / 10^6) + 1]

# ==============================================================================
# Data
# ==============================================================================
# load breakpoint data
breakpoints = fread(ARGS$breaks, sep = "\t", header = TRUE)
n_breaks = breakpoints[, .N]

# manual aggregate by locus, since you can't melt them all together at once
breakpoints_melted = rbindlist(list(
    copy(breakpoints),
    copy(breakpoints)
))
#   copy relevant locus info from second position to first
breakpoints_melted[(n_breaks + 1):(2 * n_breaks), Chrom1 := Chrom2]
breakpoints_melted[(n_breaks + 1):(2 * n_breaks), Start1 := Start2]
breakpoints_melted[(n_breaks + 1):(2 * n_breaks), End1 := End2]
breakpoints_melted[(n_breaks + 1):(2 * n_breaks), Strand1 := Strand2]
#   drop second information since it's already duplicated
breakpoints_melted[, Chrom2:= NULL]
breakpoints_melted[, Start2:= NULL]
breakpoints_melted[, End2:= NULL]
breakpoints_melted[, Strand2:= NULL]
#   fix column names
colnames(breakpoints_melted) = c(
    "Sample", "Chrom", "Start", "End",
    "Strand", "LogOdds", "Resolution"
)
# convert Chrom to ordered factor
breakpoints_melted[, Chrom := factor(Chrom, levels = CHRS, ordered = TRUE)]

# load hg38 ideogram
data(UCSC.HG38.Human.CytoBandIdeogram)

# ==============================================================================
# Analysis
# ==============================================================================
# convert breakpoints into Mbp bins and counts the rearrangements
breakpoints_melted[, StartBin := round(Start / 10^6)]
breakpoints_melted[, EndBin := round(End / 10^6)]

# convert breakpoints in genomic regions into bins
breakpoints_binned = rbindlist(lapply(
    1:breakpoints_melted[, .N],
    function(i) {
        print(i)
        start_bin = breakpoints_melted[i, StartBin]
        end_bin = breakpoints_melted[i, EndBin]
        dt = data.table(
            Sample = breakpoints_melted[i, Sample],
            Chrom = breakpoints_melted[i, Chrom],
            Bin = start_bin:end_bin,
            Count = 1
        )
        return(dt)
    }
))
# count all breakpoints in each bin
breakpoints_summed = breakpoints_binned[, sum(Count), by = c("Sample", "Chrom", "Bin")]
colnames(breakpoints_summed) = c("Sample", "Chrom", "Bin", "Count")

fwrite(
    breakpoints_summed[order(Sample, Chrom, Bin)],
    ARGS$output,
    sep = "\t",
    col.names = TRUE
)

# ==============================================================================
# Plots
# ==============================================================================
# set core component for circos plots
RCircos.Set.Core.Components(
    UCSC.HG38.Human.CytoBandIdeogram,
    NULL,
    tracks.inside = 0,
    tracks.outside = 0
)
# plot circos plots for each patient
for (s in unique(breakpoints$Sample)) {
    cat(s, "\n")
    png(
        paste0("Plots/", s, ".circos.png"),
        width = 12,
        height = 12,
        units = "cm",
        res = 300
    )
    RCircos.Set.Plot.Area()
    RCircos.Chromosome.Ideogram.Plot()
    RCircos.Link.Plot(
        link.data = breakpoints[
            Sample == s,
            .(Chrom1, Start1, End1, Chrom2, Start2, End2)
        ],
        track.num = 1,
        by.chromosome = FALSE
    )
    dev.off()
}

# number of SVs detected per patient
gg = (
    ggplot(data = breakpoints[, .N, by = Sample])
    + geom_col(aes(x = Sample, y = N, fill = Sample))
    + labs(x = "Patient", y = "Number of SVs")
    + guides(fill = FALSE)
    + theme_minimal()
    + theme(
        axis.text.x = element_text(angle = 90, vjust = 0, hjust = 0.5)
    )
)
ggsave(
    "Plots/sv-counts.png",
    height = 12,
    width = 20,
    units = "cm"
)

# number and location of SVs across all patients
gg = (
    ggplot()
    + geom_col(
        data = breakpoints_summed,
        mapping = aes(x = Bin, y = Count, fill = Sample),
        position = "stack"
    )
    # add geom_blank to ensure each facet is the correct size
    + geom_blank(data = hg38, mapping = aes(x = Bins, y = 0))
    + labs(x = "Position", y = "Count")
    + facet_grid(. ~ Chrom, scales = "free_x", space = "free_x", switch = "x")
    + theme_classic()
    + theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing = unit(0, "lines"),
        strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5)
    )
)
ggsave(
    "Plots/sv-loci.png",
    height = 12,
    width = 40,
    units = "cm"
)
