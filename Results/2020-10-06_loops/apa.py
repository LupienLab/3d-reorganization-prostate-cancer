# ==============================================================================
# Meta
# ==============================================================================
# apa.py
# --------------------------------------
# Description: Aggregate peak analysis for differential tumour-vs-benign loop calls
# Author: James Hawley

import argparse
import os.path as path
import numpy as np
import datatable as dt
from datatable import f, by, update
from cooler import Cooler
import matplotlib.pyplot as plt
import scipy.stats as stats
from tqdm import tqdm
import negspy.coordinates as nc
import logging

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

# +/- number of bins for aggregate peak analysis
SHIFT_IDX = int(RES / 1000)

# genome coordinates
hg38 = nc.get_chrominfo("hg38")

# ==============================================================================
# Functions
# ==============================================================================

# ==============================================================================
# Data
# ==============================================================================
logging.info("Loading data")
# load sample metadata
metadata = dt.fread("config.tsv", sep="\t")
# remove samples that shouldn't be included in analyses
metadata = metadata[f.Include == "Yes", :]
# extract sample IDs
SAMPLES = metadata["SampleID"].to_list()[0]
TUMOUR_SAMPLES = metadata[f.Type == "Malignant", f.SampleID].to_list()[0]
BENIGN_SAMPLES = metadata[f.Type == "Benign", f.SampleID].to_list()[0]

# load loop calls
loops = {
    "shared": dt.fread(path.join(DIR["loops"], "shared-loops.tsv"), sep="\t"),
    "tumour": dt.fread(path.join(DIR["loops"], "tumour-specific-loops.tsv"), sep="\t"),
    "benign": dt.fread(path.join(DIR["loops"], "benign-specific-loops.tsv"), sep="\t"),
}

# contact matrices
mtx = {
    s: Cooler(path.join(DIR["contacts"], s + ".mcool::/resolutions/" + str(RES)))
    for s in SAMPLES
}

# ==============================================================================
# Analysis
# ==============================================================================
# combine coordinates into UCSC-like strings for getting region extents
logging.info("Calculating loop matrix indices")
for (t, l) in loops.items():
    for i in range(l.nrows):
        # find matrix row/col IDs
        x_reg = (l[i, "chr_x"], l[i, "start_x"], l[i, "end_x"])
        y_reg = (l[i, "chr_y"], l[i, "start_y"], l[i, "end_y"])
        x_centre = int(
            sum(
                [
                    nc.chr_pos_to_genome_pos(x_reg[0], x_reg[1], hg38),
                    nc.chr_pos_to_genome_pos(x_reg[0], x_reg[2], hg38),
                ]
            )
            / (2 * RES)
        )
        y_centre = int(
            sum(
                [
                    nc.chr_pos_to_genome_pos(y_reg[0], y_reg[1], hg38),
                    nc.chr_pos_to_genome_pos(y_reg[0], y_reg[2], hg38),
                ]
            )
            / (2 * RES)
        )
        l[i, "x_left"] = x_centre - SHIFT_IDX - 1
        l[i, "x_right"] = x_centre + SHIFT_IDX
        l[i, "y_left"] = y_centre - SHIFT_IDX - 1
        l[i, "y_right"] = y_centre + SHIFT_IDX

logging.info("Aggregating matrices")
for (t, l) in loops.items():
    # set placeholder matrices
    local_loop_matrices = {
        "tumour": [np.zeros((2 * SHIFT_IDX + 1, 2 * SHIFT_IDX + 1))] * l.nrows,
        "benign": [np.zeros((2 * SHIFT_IDX + 1, 2 * SHIFT_IDX + 1))] * l.nrows,
    }
    # store local matrices in the regions of interest
    for i in tqdm(range(l.nrows)):
        local_loop_matrices["tumour"][i] = np.nanmean(
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

