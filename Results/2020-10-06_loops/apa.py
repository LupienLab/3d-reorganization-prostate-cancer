# ==============================================================================
# Meta
# ==============================================================================
# apa.py
# --------------------------------------
# Description: Aggregate peak analysis for differential tumour-vs-benign loop calls
# Author: James Hawley

import os.path as path
import numpy as np
from cooler import Cooler
import scipy.stats as stats
from tqdm import tqdm
import negspy.coordinates as nc
import logging
from cooltools import snipping
import pandas as pd
import bioframe

logging.getLogger().setLevel(logging.INFO)

# ==============================================================================
# Constants
# ==============================================================================
DIR = {
    "loops": "Loops",
    "contacts": path.join(
        "..", "..", "Data", "Processed", "2019-06-18_PCa-LowC-sequencing", "Contacts"
    ),
}

# resolution for contact matrices
RES = 5000

# +/- number of bps for aggregate peak analysis
SHIFT_SIZE = 25000

# genome coordinates
hg38 = nc.get_chrominfo("hg38")

CHROM_SIZES = bioframe.fetch_chromsizes("hg38")
CHRS = list(CHROM_SIZES.index)

# ==============================================================================
# Data
# ==============================================================================
logging.info("Loading data")
# load sample metadata
metadata = pd.read_csv("config.tsv", sep="\t")
# remove samples that shouldn't be included in analyses
metadata = metadata.loc[metadata.Include == "Yes", :]
# extract sample IDs
SAMPLES = metadata["SampleID"].tolist()
TUMOUR_SAMPLES = metadata.loc[metadata.Type == "Malignant", "SampleID"].tolist()
BENIGN_SAMPLES = metadata.loc[metadata.Type == "Benign", "SampleID"].tolist()

# load loop calls
loops = {
    "shared": pd.read_csv(path.join(DIR["loops"], "shared-loops.tsv"), sep="\t"),
    "tumour": pd.read_csv(
        path.join(DIR["loops"], "tumour-specific-loops.tsv"), sep="\t"
    ),
    "benign": pd.read_csv(
        path.join(DIR["loops"], "benign-specific-loops.tsv"), sep="\t"
    ),
}

# contact matrices
mtx = {
    s: Cooler(path.join(DIR["contacts"], s + ".mcool::/resolutions/" + str(RES)))
    for s in SAMPLES
}

# expected matrix values
expected_mtx = {
    s: pd.read_csv(path.join(DIR["contacts"], s + ".res_10000bp.exp.cis.tsv"), sep="\t")
    for s in SAMPLES
}

# ==============================================================================
# Analysis
# ==============================================================================
# combine coordinates into UCSC-like strings for getting region extents
logging.info("Calculating loop matrix indices")

# chromosomes and their sizes
supports = [(chrom, 0, CHROM_SIZES[chrom]) for chrom in CHRS]

# placeholder for all loop positions and matrix indices
windows = {
    "shared": pd.DataFrame(),
    "tumour": pd.DataFrame(),
    "benign": pd.DataFrame(),
}

for (loop_type, typed_loops) in loops.items():
    # use snipping module to get chromosome locations around each loop anchor
    windows_x = snipping.make_bin_aligned_windows(
        binsize=RES,
        chroms=typed_loops["chr_x"],
        centers_bp=(typed_loops["start_x"] + typed_loops["end_x"])
        // 2,  # integer division
        flank_bp=SHIFT_SIZE,
    )
    windows_y = snipping.make_bin_aligned_windows(
        binsize=RES,
        chroms=typed_loops["chr_y"],
        centers_bp=(typed_loops["start_y"] + typed_loops["end_y"])
        // 2,  # integer division
        flank_bp=SHIFT_SIZE,
    )
    # merge the two pairs of windows together
    windows_merged = pd.merge(
        left=windows_x,
        right=windows_y,
        left_index=True,
        right_index=True,
        suffixes=(
            "1",
            "2",
        ),  # these suffixes are required, assign_regions doesn't work properly if not '1' and '2'
    )
    # assign chromosome regions to the low/high indices
    windows[loop_type] = snipping.assign_regions(windows_merged, supports)
    windows[loop_type] = windows[loop_type].dropna()

logging.info("Aggregating matrices")
# calculate obs/exp matrix for each sample
snipper = {
    s: snipping.ObsExpSnipper(mtx[s], expected_mtx[s]) for s in SAMPLES
}
for (loop_type, typed_loops) in loops.items():
    # taken from https://cooltools.readthedocs.io/en/latest/notebooks/06_snipping-pileups.html
    stack = {
        s: snipping.pileup(windows[loop_type], snipper[s].select, snipper[s].snip) for s in SAMPLES
    }
    print(stack)
    piles = np.nanmean(stack, axis=2)
    print(piles)
    input()

    # store local matrices in the regions of interest
    for i in tqdm(range(l.nrows)):
        local_loop_matrices["tumour"][i] = np.nanconda actmean(
            [
                mtx[s].matrix()[
                    l[i, "x_left"] : l[i, "x_right"], l[i, "y_left"] : l[i, "y_right"]
                ]
                for s in TUMOUR_SAMPLES
            ],
            axis=0,
        )
        local_loop_matrices["benign"][i] = np.nanmean(
            [
                mtx[s].matrix()[
                    l[i, "x_left"] : l[i, "x_right"], l[i, "y_left"] : l[i, "y_right"]
                ]
                for s in BENIGN_SAMPLES
            ],
            axis=0,
        )
    # aggregate matrices by summing over all loops
    agg_loop_matrices = {
        sample_type: np.nansum(local_loop_matrices[sample_type], axis=0)
        for sample_type in ["tumour", "benign"]
    }
    # save data for future use
    for sample_type in ["tumour", "benign"]:
        # save list of matrices at each loop call
        np.savez_compressed(
            "Loops/{loop_type}-loops.{sample_type}-samples.npz".format(
                loop_type=t, sample_type=sample_type
            ),
            args=local_loop_matrices[sample_type],
        )
        # save aggregated matrices
        np.savez_compressed(
            "Loops/{loop_type}-loops.{sample_type}-samples.agg.npz".format(
                loop_type=t, sample_type=sample_type
            ),
            args=agg_loop_matrices[sample_type],
        )
    print(
        dt.Frame(
            {
                "Sample Type": ["Tumour", "Benign"],
                "Min": [
                    agg_loop_matrices["tumour"].min(),
                    agg_loop_matrices["benign"].min(),
                ],
                "Max": [
                    agg_loop_matrices["tumour"].max(),
                    agg_loop_matrices["benign"].max(),
                ],
            },
        )
    )

