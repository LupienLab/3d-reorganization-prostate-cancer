"""
calc-tad-similarity
==========

Calculate the BPscore for a set of samples
"""

from __future__ import division, absolute_import, print_function
import argparse
import os.path as path
import numpy as np
import pandas as pd
from tqdm import tqdm
from bpscore import bpscore

# ==============================================================================
# Constants
# ==============================================================================
TAD_DIR = path.join("resolved-TADs", "separated-TADs")
WINDOWS = list(range(3, 41))
# exclude chrY from comparisons, to account for female-derived cell lines
CHROMS = ["chr" + str(i) for i in list(range(1, 23)) + ["X"]]

# ==============================================================================
# Data
# ==============================================================================
TUMOUR_CONFIG = pd.read_csv(
    path.join("..", "..", "Data", "External", "LowC_Samples_Data_Available.tsv"),
    sep="\t"
)
TUMOUR_CONFIG = TUMOUR_CONFIG.loc[TUMOUR_CONFIG.Include == "Yes", :]

BENIGN_CONFIG = pd.read_csv(
    path.join("..", "..", "Data", "Raw", "191220_A00827_0104_AHMW25DMXX", "config.tsv"),
    sep="\t"
)
BENIGN_CONFIG = BENIGN_CONFIG.loc[BENIGN_CONFIG.Include == "Yes", :]

CELL_LINE_CONFIG = pd.read_csv(
    path.join("..", "..", "Data", "External", "Rhie_2019", "config.tsv"),
    sep="\t"
)

TUMOUR_SAMPLES = ["PCa" + str(i) for i in TUMOUR_CONFIG.loc[TUMOUR_CONFIG.Include == "Yes", "Sample ID"]]
BENIGN_SAMPLES = BENIGN_CONFIG["Sample"].tolist()
CELL_LINE_SAMPLES = CELL_LINE_CONFIG["Run_Accession"].tolist()
ALL_SAMPLES = TUMOUR_SAMPLES + BENIGN_SAMPLES + CELL_LINE_SAMPLES

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
# drop chromosome M, since this isn't included in the TAD calls
chrom_sizes.drop(index=chrom_sizes.loc[chrom_sizes.chr == "chrM", :].index, inplace=True)
# calculate total genome size
genome_size = chrom_sizes.Length.sum()

# table to store BPscore distances between sample pairs at a given window size
dists_df = pd.DataFrame(columns=["w", "s1", "s2", "dist"])

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
