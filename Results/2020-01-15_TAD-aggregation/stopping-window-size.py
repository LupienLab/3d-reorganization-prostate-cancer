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

# drop chromosome Y, since some cell lines aren't derived from males
CHROMS = ["chr" + str(i) for i in list(range(1, 23)) + ["X"]]

DIFF_THRESH = 5e-3
N_DIFFS = 1

# =============================================================================================================================
# Data
# =============================================================================================================================
print("Loading data")
metadata = pd.read_csv("config.tsv", sep="\t")
ALL_SAMPLES = metadata.loc[metadata.Include == "Yes", "SampleID"].tolist()

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
    for s in ALL_SAMPLES
}

# load chromosomes sizes
chrom_sizes = pd.read_csv(
    path.join("..", "..", "Data", "Processed", "2019-06-18_PCa-LowC-sequencing", "hg38.sizes.txt"),
    sep="\t",
    header=None,
    names=["chr", "Length"]
)
# drop chromosome Y, since some cell lines aren't derived from males
# drop chromosome M, since this isn't included in the TAD calls
chrom_sizes.drop(index=chrom_sizes.loc[chrom_sizes.chr.isin(["chrY", "chrM"]), :].index, inplace=True)
# calculate total genome size
genome_size = chrom_sizes.Length.sum()


# =============================================================================================================================
# Analysis
# =============================================================================================================================
# data to store results over window sizes
window_diffs = pd.DataFrame({
    "SampleID": ALL_SAMPLES * (len(WINDOWS) - 1),
    "w": sorted(WINDOWS[1:] * len(ALL_SAMPLES)),
    "BPscore": [0] * len(ALL_SAMPLES) * (len(WINDOWS) - 1),
    "window_step_abs_delta": [0] * len(ALL_SAMPLES) * (len(WINDOWS) - 1),
})

sup_windows = pd.DataFrame({
    "SampleID": ALL_SAMPLES,
    "w" : [0] * len(ALL_SAMPLES)
})

for s in tqdm(ALL_SAMPLES):
    # bp-score from previous window size
    b_prev = 1
    # count how many times in a row the difference threshold is met
    consistent_counter = 0
    sup_w = 0
    for w in WINDOWS[1:]:
        b = 0
        for c in CHROMS:
            # get TADs at window size w and w - 1 for chromosome c
            tads_w = tads[s][w].loc[tads[s][w].chr == c, :]
            tads_w_1 = tads[s][w - 1].loc[tads[s][w - 1].chr == c, :]
            # compare these TAD calls
            distance = bpscore(tads_w, tads_w_1)
            # scale similarity based on chromosome size
            addition = chrom_sizes.loc[chrom_sizes.chr == c, "Length"].item() / genome_size * distance
            b += addition
        # compare total BPscore (b) to BPscore for the previous window size (b_prev)
        delta = abs(b - b_prev)
        # save results
        window_diffs.loc[(window_diffs["SampleID"] == s) & (window_diffs["w"] == w), "BPscore"] = b
        window_diffs.loc[(window_diffs["SampleID"] == s) & (window_diffs["w"] == w), "window_step_abs_delta"] = delta
        # if the gain in similarity over consecutive window size steps (BPscore delta) is small enough
        if delta < DIFF_THRESH:
            # if it's not at the consistent threshold
            if consistent_counter < N_DIFFS:
                consistent_counter += 1
            # if the BPscore delta has been small enough for N_DIFFS consecutive times
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

