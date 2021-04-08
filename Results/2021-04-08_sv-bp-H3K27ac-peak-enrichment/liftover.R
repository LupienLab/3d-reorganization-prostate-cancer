# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("rtracklayer"))
suppressMessages(library("GenomicRanges"))

# ==============================================================================
# Data
# ==============================================================================
# read input file
# don't add 1 to start coords, since these come from GrCh37 VCFs
input <- fread(
    file.path("..", "..", "Data", "External", "CPC-GENE", "structural_somatic_mutation.PRAD-CA.tsv"),
    header = TRUE,
    sep = "\t"
)

sv_bkpts <- rbindlist(list(
    input[, .(
        chr = paste0("chr", chr_from),
        pos = chr_from_bkpt,
        strand = chr_from_strand
    )],
    input[, .(
        chr = paste0("chr", chr_to),
        pos = chr_to_bkpt,
        strand = chr_to_strand
    )]
))
sv_bkpts <- sv_bkpts[, .(
    chr,
    start = pos,
    end = pos,
    strand = ifelse(strand == 1, "+", "-")
)]

# convert to GRanges
gr <- GRanges(sv_bkpts)

# read liftover chain
chain <- import.chain(
    file.path("..", "..", "Data", "External", "hg19ToHg38.over.chain")
)

# ==============================================================================
# Analysis
# ==============================================================================
# lift over to hg38
gr_hg38 <- unlist(liftOver(gr, chain))

# ==============================================================================
# Save data
# ==============================================================================
gr_hg38_dt <- data.table(
    chr = as.character(seqnames(gr_hg38)),
    start = start(gr_hg38) - 1,
    end = start(gr_hg38)
)

fwrite(
    gr_hg38_dt,
    file.path("enrichment", "sv-breakpoints.bed"),
    sep = "\t",
    col.names = FALSE
)
