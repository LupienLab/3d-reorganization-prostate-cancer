# ==============================================================================
# Environment
# ==============================================================================
import numpy as np
import pandas as pd
import bioframe
import cooltools
import cooler

# ==============================================================================
# Data
# ==============================================================================
# reference genome
hg38 = bioframe.fetch_chromsizes("hg38")

# compartment calling resolution
res = 40000

# create bins of the reference genome
bins = cooler.binnify(hg38, res)

# load hg38 reference sequence
fasta_records = bioframe.load_fasta("data/hg38.fa")

# ==============================================================================
# Analysis
# ==============================================================================
# calculate GC content in each bin
bins["GC"] = bioframe.tools.frac_gc(bins, fasta_records)

# ==============================================================================
# Save data
# ==============================================================================
bins.to_csv(
	"gc-content-phase.bedGraph",
	sep="\t",
	names=True,
)
