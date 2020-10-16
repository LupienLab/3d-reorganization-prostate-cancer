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
import pickle

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

CHROM_SIZES = pd.read_csv(
    path.join(DIR["contacts"], "..", "hg38.sizes.txt"),
    sep="\t",
    index_col=0,
    header=None,
    names=["size"],
)
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
supports = [(chrom, 0, CHROM_SIZES.at[chrom, "size"]) for chrom in CHRS]

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
# save serialized object
snipper_obj = open("Loops/snipper.obj", "wb")
pickle.dump(snipper, snipper_obj)
snipper_obj.close()

for (loop_type, typed_loops) in tqdm(loops.items()):
    logging.info("Loops {}".format(loop_type))
    # taken from https://cooltools.readthedocs.io/en/latest/notebooks/06_snipping-pileups.html
    # create a stack of obs/exp matrices based on the locations in windows[loop_type]
    print(windows[loop_type])
    stack = {
        s: snipping.pileup(windows[loop_type], snipper[s].select, snipper[s].snip) for s in SAMPLES
    }
    # save serialized object
    stack_obj = open(".".join(["Loops/stack", loop_type, "obj"]), "wb")
    pickle.dump(stack, stack_obj)
    stack_obj.close()
    # sum each stack to create an obs/exp pileup for each sample
    piles = {s: np.nanmean(stack[s], axis=2) for s in SAMPLES}
    # save serialized object
    piles_obj = open(".".join(["Loops/pile", loop_type, "obj"]), "wb")
    pickle.dump(piles, piles_obj)
    piles_obj.close()
    # take the mean over each condition (benign/tumour sample)
    condition_piles = {
        "tumour": np.nanmean([piles[s] for s in TUMOUR_SAMPLES], axis=0),
        "benign": np.nanmean([piles[s] for s in BENIGN_SAMPLES], axis=0),
    }
    # save serialized object
    condition_piles_obj = open(".".join(["Loops/condition_pile", loop_type, "obj"]), "wb")
    pickle.dump(condition_piles, condition_piles_obj)
    condition_piles_obj.close()
    print(
        pd.DataFrame(
            {
                "Sample Type": ["Tumour", "Benign"],
                "Min": [
                    condition_piles["tumour"].min(),
                    condition_piles["benign"].min(),
                ],
                "Max": [
                    condition_piles["tumour"].max(),
                    condition_piles["benign"].max(),
                ],
            },
        )
    )


