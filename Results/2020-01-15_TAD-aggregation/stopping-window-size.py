"""
stopping-window-size
==========

Calculate the optimal window size to stop calling TADs at
"""

from __future__ import division, absolute_import, print_function
import argparse
import os.path as path
import numpy as np
import pandas as pd
from bpscore import bpscore
from tqdm import tqdm

# =============================================================================================================================
# Constants
# =============================================================================================================================
TAD_DIR = path.join("resolved-TADs", "separated-TADs")
WINDOWS = list(range(3, 41))
CHROMS = ["chr" + str(i) for i in list(range(1, 23)) + ["X", "Y"]]
DIFF_THRESH = 5e-3
N_DIFFS = 2

# =============================================================================================================================
# Data
# =============================================================================================================================
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


# =============================================================================================================================
# Analysis
# =============================================================================================================================
# data to store results over window sizes
window_diffs = pd.DataFrame({
    "SampleID": SAMPLES * (len(WINDOWS) - 1),
    "w": sorted(WINDOWS[1:] * len(SAMPLES)),
    "diff": [0] * len(SAMPLES) * (len(WINDOWS) - 1),
    "abs_delta": [0] * len(SAMPLES) * (len(WINDOWS) - 1),
})

sup_windows = pd.DataFrame({
    "SampleID": SAMPLES,
    "w" : [0] * len(SAMPLES)
})

for s in tqdm(SAMPLES):
    b_prev = 1
    consistent_counter = 0
    sup_w = 0
    for w in WINDOWS[1:]:
        b = 0
        for c in CHROMS:
            tads_w = tads[s][w].loc[tads[s][w].chr == c, :]
            tads_w_1 = tads[s][w - 1].loc[tads[s][w - 1].chr == c, :]
            distance = bpscore(tads_w, tads_w_1)
            addition = chrom_sizes.loc[chrom_sizes.chr == c, "Length"].item() / genome_size * distance
            b += addition
        delta = abs(b - b_prev)
        # save results
        window_diffs.loc[(window_diffs["SampleID"] == s) & (window_diffs["w"] == w), "diff"] = b
        window_diffs.loc[(window_diffs["SampleID"] == s) & (window_diffs["w"] == w), "abs_delta"] = delta
        if delta < DIFF_THRESH:
            if consistent_counter < N_DIFFS:
                consistent_counter += 1
            elif sup_w == 0:
                sup_w = w
                sup_windows.loc[(sup_windows["SampleID"] == s), "w"] = sup_w
        else:
            consistent_counter = 0
        b_prev = b


# =============================================================================================================================
# Save data
# =============================================================================================================================
window_diffs.to_csv(
    path.join("Statistics", "tad-similarity-deltas.tsv"),
    sep="\t",
    index=False
)
sup_windows.to_csv(
    path.join("Statistics", "tad-similarity-supremum.tsv"),
    sep="\t",
    index=False
)
