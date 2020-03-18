"""
create-indicued-regions
==========

Create induced overlaps from all sample TADs into a BED file
"""

import os.path as path
import numpy as np
import pandas as pd
from tqdm import tqdm

# ==============================================================================
# Constants
# ==============================================================================
TAD_DIR = path.join(
    "..", "2020-01-15_TAD-aggregation", "resolved-TADs"
)

# window size to check TADs for
W = 3

# chromosomes
CHRS = ["chr" + str(i) for i in list(range(1, 23)) + ["X", "Y"]]

# ==============================================================================
# Data
# ==============================================================================
# load metadata
CONFIG = pd.read_csv(
    path.join("..", "..", "Data", "External", "LowC_Samples_Data_Available.tsv"),
    sep="\t",
    index_col=False,
)
SAMPLES = ["PCa" + str(i) for i in CONFIG["Sample ID"]]

# load TADs
print("Loading TADs")
bounds = {
    s: pd.read_csv(
        path.join(TAD_DIR, s + ".40000bp.aggregated-boundaries.tsv"),
        sep="\t",
        header=0,
    )
    for s in SAMPLES
}
# only keep the w=3 boundaries
bounds = {s: bounds[s].loc[bounds[s].w.str.contains(str(W) + "\|?"), :] for s in SAMPLES}


# ==============================================================================
# Analysis
# ==============================================================================
# pre-allocate memory for TADs instead of appending each time
# (later I'll delete all the rows with None's)`
regions = pd.DataFrame({
    "chr": [None] * sum(bounds[s].shape[0] for s in SAMPLES),
    "start": [None] * sum(bounds[s].shape[0] for s in SAMPLES),
    "end": [None] * sum(bounds[s].shape[0] for s in SAMPLES)
})

n_regions = 0

# iterate over each chromosome
for c in CHRS:
    i = {s: 0 for s in SAMPLES}
    chr_bounds = {s: bounds[s].loc[bounds[s].chr == c, :] for s in SAMPLES}
    lower = 0
    # step over each domain boundary
    # while the counter is always in bounds for each sample
    pbar = tqdm(unit="bound")
    while np.all([i[s] < chr_bounds[s].shape[0] for s in SAMPLES]):
        regions.loc[n_regions, "chr"] = c
        regions.loc[n_regions, "start"] = lower
        mindex = np.argmin([chr_bounds[s].iloc[i[s]].pos for s in SAMPLES])
        min_sample = SAMPLES[mindex]
        upper = chr_bounds[min_sample].iloc[i[min_sample]].pos 
        # only record if there is distance between the lower and upper bounds
        if upper > lower:
            regions.loc[n_regions, "end"] = upper
            lower = upper
        i[min_sample] += 1
        n_regions += 1
        pbar.update()

# remove regions with None as end
regions.dropna(inplace=True)

# ==============================================================================
# Save data
# ==============================================================================
regions.to_csv("TAD-induced-regions.bed", index=False, sep="\t", header=False)
