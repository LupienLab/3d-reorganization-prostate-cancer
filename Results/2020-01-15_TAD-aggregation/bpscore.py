"""
bpscore
==========

Calculate the BPscore for a set of samples
"""

from __future__ import division, absolute_import, print_function
import argparse
import os.path as path
import numpy as np
import pandas as pd
from tqdm import tqdm

# ==============================================================================
# Constants
# ==============================================================================
TAD_DIR = path.join("resolved-TADs", "separated-TADs")
WINDOWS = list(range(3, 41))
CHROMS = ["chr" + str(i) for i in list(range(1, 23)) + ["X", "Y"]]

# ==============================================================================
# Functions
# ==============================================================================
def bpscore(tads_a, tads_b, lower=None, upper=None):
    """
    Calculate the BPscore for two sets of TADs

    Parameters
    ----------
    tads_a, tads_b :
        TADs to compare
    lower : int
        Lower bound to consider
    upper : int
        Upper bound to consider
    """
    if lower is None:
        lower = max(tads_a["start"].min(), tads_b["start"].min())
    if upper is None:
        upper = min(tads_a["end"].max(), tads_b["end"].max())
    N = upper - lower
    a = sorted(tads_a["start"].tolist() + tads_a["end"].tolist())
    # winsorize the ends so they match upper and lower
    a = np.unique([min(max(lower, x), upper) for x in a])
    b = sorted(tads_b["start"].tolist() + tads_b["end"].tolist())
    b = np.unique([min(max(lower, x), upper) for x in b])
    # counters for a, b, repsectively
    i, j = (1, 1)
    # similarity
    s = 0
    while i < len(a) and j < len(b):
        overlap = min(a[i], b[j]) - max(a[i - 1], b[j - 1])
        s += overlap ** 2 / max(a[i] - a[i - 1], b[j] - b[j - 1])
        if b[j] > a[i]:
            i += 1
        else:
            j += 1
    return 1 - s / N

# ==============================================================================
# Data
# ==============================================================================
CONFIG = pd.read_csv(
    path.join("..", "..", "Data", "External", "LowC_Samples_Data_Available.tsv"),
    sep="\t",
    index_col=False,
)
SAMPLES = ["PCa" + str(i) for i in CONFIG["Sample ID"]]

# load TADs for each patient
tads = {
    s: {
        w: pd.read_csv(
            path.join(TAD_DIR, s + ".40000bp.w_" + str(w) + ".domains.bed"),
            sep="\t",
            names=["chr", "start", "end", "lower_persistence", "upper_persistence",],
        )
        for w in WINDOWS
    }
    for s in SAMPLES
}

# load chromosomes sizes
chrom_sizes = pd.read_csv(
    path.join("..", "..", "Data", "Processed", "2019-06-18_PCa-LowC-sequencing", "hg38.sizes.txt"),
    sep="\t",
    header=None,
    names=["chr", "Length"]
)
# drop chromosome M, since this isn't included in the TAD calls
chrom_sizes.drop(index=chrom_sizes.loc[chrom_sizes.chr == "chrM", :].index, inplace=True)
# calculate total genome size
genome_size = chrom_sizes.Length.sum()

# table to store BPscore distances between sample pairs at a given window size
dists_df = pd.DataFrame(columns=["w", "s1", "s2", "dist"])

for w in tqdm(WINDOWS, desc="Window sizes", position=0):
    for i, si in tqdm(enumerate(SAMPLES), total=len(SAMPLES), desc="Sample i", position=1, leave=False):
        for j, sj in tqdm(enumerate(SAMPLES), total=len(SAMPLES), desc="Sample j", position=2, leave=False):
            if j <= i:
                continue
            # bpscore between these two samples
            b = 0
            for c in CHROMS:
                tads_i = tads[si][w].loc[tads[si][w].chr == c, :]
                tads_j = tads[sj][w].loc[tads[sj][w].chr == c, :]
                distance = bpscore(tads_i, tads_j)
                addition = chrom_sizes.loc[chrom_sizes.chr == c, "Length"].item() / genome_size * distance
                b += addition
            # store BPscore in matrix for this window size
            dists_df = dists_df.append({"w": w, "s1": si, "s2": sj, "dist": b}, ignore_index=True, sort=False)
            dists_df = dists_df.append({"w": w, "s1": sj, "s2": si, "dist": b}, ignore_index=True, sort=False)


# =============================================================================================================================
# Save data
# =============================================================================================================================
dists_df.to_csv(
    path.join("Statistics", "tad-distances.tsv"),
    sep="\t",
    index=False
)
