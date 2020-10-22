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
loops = pd.read_csv(
    path.join(DIR["loops"], "merged-loops.sample-counts.tsv"),
    sep="\t",
)

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
snipper = {s: snipping.ObsExpSnipper(mtx[s], expected_mtx[s]) for s in SAMPLES}
# save serialized object
snipper_obj = open("Loops/snipper.obj", "wb")
pickle.dump(snipper, snipper_obj)
snipper_obj.close()

logging.info("Aggregating matrices at loops")
# taken from https://cooltools.readthedocs.io/en/latest/notebooks/06_snipping-pileups.html
# create a stack of obs/exp matrices based on the locations in windows
# this is the part that takes the longest time
stack = {
    s: snipping.pileup(windows, snipper[s].select, snipper[s].snip)
    for s in SAMPLES
}
# save serialized object
stack_obj = open("Loops/stack.obj", "wb")
pickle.dump(stack, stack_obj)
stack_obj.close()

logging.info("Done")

