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

SAMPLE_TYPES = ["benign", "tumour"]

# resolution for contact matrices
RES = 10000

# +/- number of bps for aggregate peak analysis
SHIFT_SIZE = 300000

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
# Functions
# ==============================================================================
def filter_arr(x: np.ndarray) -> np.ndarray:
    return x[~np.isnan(x) & np.isfinite(x)]


# ==============================================================================
# Data
# ==============================================================================
logging.info("Loading data")
# load sample metadata
metadata = pd.read_csv("config.tsv", sep="\t")
# remove samples that shouldn't be included in analyses
metadata = metadata.loc[metadata.Include == "Yes", :]

# extract sample IDs
SAMPLES = {
    "all": metadata.loc[metadata.Include == "Yes", "SampleID"].tolist(),
    "tumour": metadata.loc[metadata.Type == "Malignant", "SampleID"].tolist(),
    "benign": metadata.loc[metadata.Type == "Benign", "SampleID"].tolist(),
}

# load loop calls
loops = pd.read_csv(
    path.join(DIR["loops"], "merged-loops.sample-counts.tsv"),
    sep="\t",
)
# restrict to loops detected in >= 2 samples
loops = loops.loc[loops["Benign"] + loops["Malignant"] >= 2, :]

# contact matrices
mtx = {
    s: Cooler(path.join(DIR["contacts"], s + ".mcool::/resolutions/" + str(RES)))
    for s in SAMPLES["all"]
}

# expected matrix values
expected_mtx = {
    s: pd.read_csv(path.join(DIR["contacts"], s + ".res_10000bp.exp.cis.tsv"), sep="\t")
    for s in SAMPLES["all"]
}

# ==============================================================================
# Analysis
# ==============================================================================
# combine coordinates into UCSC-like strings for getting region extents
logging.info("Calculating loop matrix indices")

# chromosomes and their sizes
supports = [(chrom, 0, CHROM_SIZES.at[chrom, "size"]) for chrom in CHRS]

# use snipping module to get chromosome locations around each loop anchor
windows_x = snipping.make_bin_aligned_windows(
    binsize=RES,
    chroms=loops["chr_x"],
    centers_bp=(loops["start_x"] + loops["end_x"]) // 2,  # integer division
    flank_bp=SHIFT_SIZE,
)
windows_y = snipping.make_bin_aligned_windows(
    binsize=RES,
    chroms=loops["chr_y"],
    centers_bp=(loops["start_y"] + loops["end_y"]) // 2,  # integer division
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
windows = snipping.assign_regions(windows_merged, supports)
windows = windows.dropna()

logging.info("Calculating Obs/Exp matrices")
# calculate obs/exp matrix for each sample
snipper = {s: snipping.ObsExpSnipper(mtx[s], expected_mtx[s]) for s in SAMPLES["all"]}

logging.info("Aggregating matrices at loops")
# taken from https://cooltools.readthedocs.io/en/latest/notebooks/06_snipping-pileups.html
# create a stack of obs/exp matrices based on the locations in windows
# this is the part that takes the longest time
stack = {
    s: snipping.pileup(windows, snipper[s].select, snipper[s].snip)
    for s in SAMPLES["all"]
}

logging.info("Summing over condition-specific interactions")
# sum over groups for each locus to rank the differential loop enrichment across loci
conditional_stack = {
    # 3. coerce into a 3D array: (loop locus, row, col)
    sample_type: np.array(
        [
            # 2. get the sample-type mean of the obs/exp matrices for the given locus
            np.nanmean(
                # 1. get all obs/exp matrices from samples of the specific type for a given locus, i
                [stack[s][:, :, i] for s in SAMPLES[sample_type]],
                axis=0,
            )
            for i in range(stack["PCa13266"].shape[2])
        ]
    )
    for sample_type in SAMPLE_TYPES
}


# get differential obs/exp matrices for each loop locus
# a 3D array: (loop locus, row, col)
logging.info(
    "Calculating difference between benign and tumour samples in condition-specific interactions"
)
conditional_stack_differential = np.array(
    [
        #  calculate differential between tumour and benign for a given locus, i
        np.nanmean(
            filter_arr(
                np.log2(
                    conditional_stack["tumour"][i, :, :]
                    / conditional_stack["benign"][i, :, :]
                )
            )
        )
        for i in range(stack["PCa13266"].shape[2])
    ]
)

# sum over the stack of obs/exp values over each loop call
pileup = {s: np.nanmean(stack[s], axis=2) for s in SAMPLES["all"]}

# sum over each tumour/benign sample, in groups
conditional_pileup = {
    sample_type: np.nanmean([pileup[s] for s in SAMPLES[sample_type]], axis=0)
    for sample_type in SAMPLE_TYPES
}

# differential testing between tumour and benign samples
conditional_differential = conditional_pileup["tumour"] - conditional_pileup["benign"]


# ==============================================================================
# Save data
# ==============================================================================
logging.info("Saving serialized data")

# save snipper
snipper_obj = open(path.join(DIR["loops"], "snipper.obj"), "wb")
pickle.dump(snipper, snipper_obj)
snipper_obj.close()

# save stack
stack_obj = open(path.join(DIR["loops"], "stack.obj"), "wb")
pickle.dump(stack, stack_obj)
stack_obj.close()

# save stack in each condition
cdtnl_stack_obj = open(path.join(DIR["loops"], "stack.conditional.obj"), "wb")
pickle.dump(conditional_stack, cdtnl_stack_obj)
cdtnl_stack_obj.close()

# save differential stack object
difftl_stack_obj = open(
    path.join(DIR["loops"], "stack.conditional.differential.obj"), "wb"
)
pickle.dump(conditional_stack_differential, difftl_stack_obj)
difftl_stack_obj.close()

# save pileups
pileup_obj = open(path.join(DIR["loops"], "pileup.obj"), "wb")
pickle.dump(pileup, pileup_obj)
pileup_obj.close()

# save conditional pileup
cdtnl_pileup_obj = open(path.join(DIR["loops"], "pileup.conditional.obj"), "wb")
pickle.dump(conditional_pileup, cdtnl_pileup_obj)
cdtnl_pileup_obj.close()

# save differential contacts
difftl_obj = open(path.join(DIR["loops"], "differential.obj"), "wb")
pickle.dump(conditional_differential, difftl_obj)
difftl_obj.close()
