# ==============================================================================
# Environment
# ==============================================================================
print("Loading packages")
import pandas as pd
import bioframe
import cooler

# ==============================================================================
# Data
# ==============================================================================
print("Loading reference genome")
# reference genome
hg38 = bioframe.fetch_chromsizes("hg38")

# load hg38 reference sequence
fasta_records = bioframe.load_fasta(
    "/mnt/work1/data/genomes/human/hg38/iGenomes/Sequence/WholeGenomeFasta/genome.fa"
)

print("Binning the genome")
# compartment calling resolution
res = 40000

# create bins of the reference genome
bins = cooler.binnify(hg38, res)

# ==============================================================================
# Analysis
# ==============================================================================
print("Calculating GC content")
# calculate GC content in each bin
bins["GC"] = bioframe.tools.frac_gc(bins, fasta_records)

# ==============================================================================
# Save data
# ==============================================================================
print("Saving bedGraph")
bins.to_csv(
    "gc-content-phase.bedGraph",
    sep="\t",
    index=False,
    header=True,
)
