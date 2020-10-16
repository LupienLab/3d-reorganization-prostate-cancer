# ==============================================================================
# Meta
# ==============================================================================
# plot-apa
# --------------------------------------
# Description: Plot aggregate peak analyses results
# Author: James Hawley

import logging

import os.path as path
import numpy as np
import matplotlib as mpl

mpl.use("agg")

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import scipy.stats as stats
import pickle
from datatable import dt, f, update, by

logging.getLogger().setLevel(logging.INFO)


# ==============================================================================
# Constants
# ==============================================================================
LOOP_DIR = "Loops"

# resolution for contact matrices
RES = 5000

# +/- number of bps for aggregate peak analysis
SHIFT_SIZE = 25000

LOOP_TYPES = ["shared", "benign", "tumour"]
SAMPLE_TYPES = ["benign", "tumour"]

# ==============================================================================
# Functions
# ==============================================================================

# ==============================================================================
# Data
# ==============================================================================
logging.info("Loading data")
# load metadata
metadata = dt.fread("config.tsv", sep="\t")
SAMPLES = metadata[f.Include == "Yes", f.SampleID].to_list()[0]

# load pileups
loop_matrices = {
    s: {
        loop_type: pickle.load(
            open(
                path.join(LOOP_DIR, "pile.{loop_type}.obj".format(loop_type=loop_type)),
                "rb",
            )
        )[s]
        for loop_type in LOOP_TYPES
    }
    for s in SAMPLES
}

# ==============================================================================
# Plots
# ==============================================================================
logging.info("Plotting")
# create grid for plotting
gs = GridSpec(
    nrows=len(LOOP_TYPES),
    ncols=len(SAMPLES) + 1,
    width_ratios=[20] * len(SAMPLES) + [1],
)

# plotting options
plt.figure(figsize=(4 * len(SAMPLES), 20))
opts = dict(
    # vmin=-0.75,
    # vmax=0.75,
    extent=[-10, 10, -10, 10],
    cmap="coolwarm",
)

# create heatmap for each plot
for i, loop_type in enumerate(LOOP_TYPES):
    for j, s in enumerate(SAMPLES):
        ax = plt.subplot(gs[i, j])
        img = ax.matshow(np.log2(loop_matrices[s][loop_type]), **opts)
        ax.xaxis.tick_bottom()
        ax.set_xlabel(s.title())
        ax.set_ylabel(loop_type.title())
        if i < 2:
            ax.xaxis.set_visible(False)
        if j > 0:
            ax.yaxis.set_visible(False)

ax = plt.subplot(gs[:, len(SAMPLES)])
ax.set_xlabel("log2(Obs / Exp)\nContact Frequency")
ax.yaxis.tick_right()
plt.colorbar(img, cax=ax)
plt.savefig("Plots/apa.by-sample.png")
