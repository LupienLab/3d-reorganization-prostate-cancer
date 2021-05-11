"""
stopping-window-size
==========

Calculate the optimal window size to stop calling TADs at
"""

import argparse
import os.path as path
import numpy as np
import pandas as pd
from bpscore import bpscore
from tqdm import tqdm

PARSER = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)
PARSER.add_argument(
    "-c", "--count",
    type=int,
    help="Downsampling count number",
    default=300000000,
)
PARSER.add_argument(
    "-r", "--res",
    type=int,
    help="Contact matrix resolution",
    default=40000,
)
PARSER.add_argument(
    "-i", "--in-dir",
    type=str,
    help="Input directory containing TAD calls",
    default="TADs",
)
PARSER.add_argument(
    "-s", "--sample-id",
    type=str,
    help="Sample ID to include",
    nargs="+",
    dest="sample_ids",
    required=True,
)
ARGS = PARSER.parse_args()

# =============================================================================================================================
# Constants
# =============================================================================================================================
TAD_DIR = ARGS.in_dir
WINDOWS = list(range(2, 41))
RESOLUTION = ARGS.res
DWNSAMPLE = ARGS.count
ALL_SAMPLES = ARGS.sample_ids

# drop chromosome Y, since some cell lines aren't derived from males
CHROMS = ["chr" + str(i) for i in list(range(1, 23)) + ["X"]]

DIFF_THRESH_1o = 8e-2
DIFF_THRESH_2o = 1e-2
N_DIFFS = 2

# =============================================================================================================================
# Data
# =============================================================================================================================
print("Loading data")

# load TADs for each patient
tads = {
    s: {
        w: pd.read_csv(
            path.join(TAD_DIR, "{s}.{c}.res_{r}bp.window_{w}.domains.bed".format(s=s, c=DWNSAMPLE, r=RESOLUTION, w=w)),
            sep="\t",
            names=["chr", "start", "end", "type"],
        ) for w in WINDOWS
    } for s in ALL_SAMPLES
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

# bp-score from previous window size
b_prev = {s: 1 for s in ALL_SAMPLES}
# count how many times in a row the difference threshold is met
consistent_counter = 0
sup_w = 0
for w in tqdm(WINDOWS[1:]):
    # counter for BPscore
    b = {s: 0 for s in ALL_SAMPLES}
    for s in ALL_SAMPLES:
        for c in CHROMS:
            # get TADs at window size w and w - 1 for chromosome c
            tads_w = tads[s][w].loc[tads[s][w].chr == c, :]
            tads_w_1 = tads[s][w - 1].loc[tads[s][w - 1].chr == c, :]
            # compare these TAD calls
            distance = bpscore(tads_w, tads_w_1)
            # scale similarity based on chromosome size
            addition = chrom_sizes.loc[chrom_sizes.chr == c, "Length"].item() / genome_size * distance
            b[s] += addition
        # compare total BPscore (b) to BPscore for the previous window size (b_prev)
        delta = abs(b[s] - b_prev[s])
        # save results
        window_diffs.loc[(window_diffs["SampleID"] == s) & (window_diffs["w"] == w), "BPscore"] = b[s]
        window_diffs.loc[(window_diffs["SampleID"] == s) & (window_diffs["w"] == w), "window_step_abs_delta"] = delta
#     print(b)
#     print(b_prev)
#     print(window_diffs.loc[window_diffs.w == w, :])
#     input()
    b_prev = b
    # if the similarity between window sizes is small
    print(window_diffs.loc[window_diffs.w == w, "BPscore"].max(), window_diffs.loc[window_diffs.w == w, "window_step_abs_delta"].max())
    if window_diffs.loc[window_diffs.w == w, "BPscore"].max() <= DIFF_THRESH_1o:
        # if the gain in similarity over consecutive window size steps is small enough
        if window_diffs.loc[window_diffs.w == w, "window_step_abs_delta"].max() <= DIFF_THRESH_2o:
            # if this is the n-th time in a row that this has happened
            if consistent_counter == N_DIFFS:
                # if this is the first time this similarity threshold has been reached
                if sup_w == 0:
                    sup_w = w
                    sup_windows.loc[:, "w"] = sup_w
            else:
                consistent_counter += 1
        else:
            consistent_counter = 0
    else:
        consistent_counter = 0

print(sup_w)


# =============================================================================================================================
# Save data
# =============================================================================================================================
window_diffs.to_csv(
    path.join("Statistics", "tad-similarity-deltas." + str(DWNSAMPLE) + ".tsv"),
    sep="\t",
    index=False
)
sup_windows.to_csv(
    path.join("Statistics", "tad-similarity-supremum." + str(DWNSAMPLE) + ".tsv"),
    sep="\t",
    index=False
)

