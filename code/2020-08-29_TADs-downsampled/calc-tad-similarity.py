"""
calc-tad-similarity
==========

Calculate the BPscore for a set of samples
"""

import argparse
import os.path as path
import numpy as np
import pandas as pd
from tqdm import tqdm
from bpscore import bpscore

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
    "-m", "--min",
    type=int,
    help="Minimum window size to consider",
    default=3,
)
PARSER.add_argument(
    "-M", "--max",
    type=int,
    help="Maximum window size to consider",
    default=20,
)
PARSER.add_argument(
    "-s", "--sample-id",
    type=str,
    help="Sample ID(s) to include",
    nargs="+",
    dest="sample_ids",
    required=True,
)
ARGS = PARSER.parse_args()


# ==============================================================================
# Constants
# ==============================================================================
TAD_DIR = path.join("Aggregated-TADs", "separated-TADs")
RESOLUTION = 40000
WINDOWS = list(range(ARGS.min, ARGS.max + 1))
DWNSAMPLE = ARGS.count
ALL_SAMPLES = ARGS.sample_ids

# exclude chrY from comparisons, to account for female-derived cell lines
CHROMS = ["chr" + str(i) for i in list(range(1, 23)) + ["X"]]

# ==============================================================================
# Data
# ==============================================================================
# load TADs for each patient
tads = {
    s: {
        w: pd.read_csv(
            path.join(TAD_DIR, s + "." + str(DWNSAMPLE) + ".res_" + str(RESOLUTION) + "bp.window_" + str(w) + ".domains.tsv"),
            sep="\t",
            header=None,
            names=["chr", "start", "end", "persistence_left", "persistence_right", "type"],
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
# drop chromosome M, since this isn't included in the TAD calls
chrom_sizes.drop(index=chrom_sizes.loc[chrom_sizes.chr == "chrM", :].index, inplace=True)
# calculate total genome size
genome_size = chrom_sizes.Length.sum()

# table to store BPscore distances between sample pairs at a given window size
dists_df = pd.DataFrame(columns=["w", "s1", "s2", "BPscore"])

for w in tqdm(WINDOWS, desc="Window sizes", position=0):
    for i, si in tqdm(enumerate(ALL_SAMPLES), total=len(ALL_SAMPLES), desc="Sample i", position=1, leave=False):
        for j, sj in tqdm(enumerate(ALL_SAMPLES), total=len(ALL_SAMPLES), desc="Sample j", position=2, leave=False):
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
            dists_df = dists_df.append({"w": w, "s1": si, "s2": sj, "BPscore": b}, ignore_index=True, sort=False)
            dists_df = dists_df.append({"w": w, "s1": sj, "s2": si, "BPscore": b}, ignore_index=True, sort=False)


# =============================================================================================================================
# Save data
# =============================================================================================================================
dists_df.to_csv(
    path.join("Statistics", "tad-distances." + str(DWNSAMPLE) + ".tsv"),
    sep="\t",
    index=False
)
