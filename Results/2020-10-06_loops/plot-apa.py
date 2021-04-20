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
import pandas as pd
from itertools import chain

logging.getLogger().setLevel(logging.INFO)


# ==============================================================================
# Constants
# ==============================================================================
LOOP_DIR = "Loops"

# resolution for contact matrices
RES = 10000

# +/- number of bps for aggregate peak analysis
SHIFT_SIZE = 300000

SAMPLE_TYPES = ["benign", "tumour"]


# ==============================================================================
# Functions
# ==============================================================================
def filter_arr(x: np.ndarray) -> np.ndarray:
    return x[~np.isnan(x) & np.isfinite(x)]


# ==============================================================================
# Data
# ==============================================================================
logging.info("Loading data")
# load metadata
metadata = pd.read_csv("config.tsv", sep="\t")
SAMPLES = {
    "all": metadata.loc[metadata.Include == "Yes", "SampleID"].tolist(),
    "tumour": metadata.loc[metadata.Type == "Malignant", "SampleID"].tolist(),
    "benign": metadata.loc[metadata.Type == "Benign", "SampleID"].tolist(),
}

# load loop calls
loops = pd.read_csv("Loops/merged-loops.sample-counts.tsv", sep="\t")
tumour_loop_idx = loops[loops.Benign == 0].index
benign_loop_idx = loops[loops.Malignant == 0].index
shared_loop_idx = loops[(loops.Malignant > 0) & (loops.Benign > 0)].index

pileup = pickle.load(open(path.join(LOOP_DIR, "pileup.obj"), "rb"))
conditional_stack = pickle.load(
    open(path.join(LOOP_DIR, "stack.conditional.obj"), "rb")
)
conditional_pileup = pickle.load(
    open(path.join(LOOP_DIR, "pileup.conditional.obj"), "rb")
)
conditional_differential = pickle.load(
    open(path.join(LOOP_DIR, "differential.obj"), "rb")
)


# ==============================================================================
# Analysis
# ==============================================================================
# rank the loops by their mean differential between benign and tumour samples
difference_stack = np.array(
    [
        conditional_stack["tumour"][i, :, :] - conditional_stack["benign"][i, :, :]
        for i in range(conditional_stack["benign"].shape[0])
    ]
)

normed_stack = np.array(
    [
        np.linalg.norm(
            conditional_stack["tumour"][i, :, :] - conditional_stack["benign"][i, :, :],
            ord="fro",  # Frobenius norm on the difference matrix
        )
        for i in range(conditional_stack["benign"].shape[0])
    ]
)

# replace any NaN's with 0 to not mess with ranking
normed_stack[np.isnan(normed_stack)] = 0
# rank the differential so we can identify the strongest "benign-specific" and "tumour-specific" loops
loop_ordering = np.argsort(normed_stack)

n_top_loops = 10
top_loops_idx = list(
    chain(
        loop_ordering[0 : (n_top_loops // 2)],
        loop_ordering[-(n_top_loops // 2) :],
    )
)


# calculate APA for "tumour-specific" and "benign-specific" loops
specific_apa = {
    "tumour": {
        "tumour-specific": np.nanmean(
            conditional_stack["tumour"][tumour_loop_idx, :, :], axis=0
        ),
        "benign-specific": np.nanmean(
            conditional_stack["tumour"][benign_loop_idx, :, :], axis=0
        ),
        "shared": np.nanmean(
            conditional_stack["tumour"][shared_loop_idx, :, :], axis=0
        ),
    },
    "benign": {
        "tumour-specific": np.nanmean(
            conditional_stack["benign"][tumour_loop_idx, :, :], axis=0
        ),
        "benign-specific": np.nanmean(
            conditional_stack["benign"][benign_loop_idx, :, :], axis=0
        ),
        "shared": np.nanmean(
            conditional_stack["benign"][shared_loop_idx, :, :], axis=0
        ),
    },
}

# ==============================================================================
# Plots
# ==============================================================================
# create heatmap for each pileup plot
ncols = len(SAMPLES["all"]) + 1
nrows = 1
# create grid specification
gs = GridSpec(
    nrows=nrows,
    ncols=ncols,
    width_ratios=[20] * (ncols - 1) + [1],
)
# plotting options
plt.figure(figsize=(4 * (ncols - 1), 4 * nrows))
opts = dict(
    extent=[-10, 10, -10, 10],
    cmap="coolwarm",
)  # vmin=-2, vmax=2,)
# make component plots
for j, s in enumerate(SAMPLES["all"]):
    ax = plt.subplot(gs[0, j])
    img = ax.matshow(np.log2(pileup[s]), **opts)
    # add x axis labels to bottom-most subplots
    ax.set_xlabel(s)
    ax.xaxis.set_visible(True)
    ax.xaxis.tick_bottom()


# add colourbar
ax = plt.subplot(gs[:, ncols - 1])
ax.set_xlabel("log2(Obs / Exp)\nContact Frequency")
ax.yaxis.tick_right()
plt.colorbar(img, cax=ax)
plt.savefig("Plots/apa.pileup.samples.png")
plt.close()


# create heatmap for each conditional pileup plot
# get rows and columns per loop type
ncols = len(SAMPLE_TYPES) + 1
nrows = 1
# create grid specification
gs = GridSpec(
    nrows=nrows,
    ncols=ncols,
    width_ratios=[20] * (ncols - 1) + [1],
)
# plotting options
plt.figure(figsize=(4 * (ncols - 1), 4 * nrows))
# make component plots
opts = dict(
    extent=[
        -(SHIFT_SIZE // 1000),  # values in kbp
        (SHIFT_SIZE // 1000),
        -(SHIFT_SIZE // 1000),
        (SHIFT_SIZE // 1000),
    ],
    cmap="coolwarm",
)
for j, sample_type in enumerate(SAMPLE_TYPES):
    ax = plt.subplot(gs[0, j])
    img = ax.matshow(np.log2(conditional_pileup[sample_type]), **opts)
    # add x axis labels to bottom-most subplots
    ax.set_xlabel(sample_type.title())
    ax.xaxis.set_visible(True)
    ax.xaxis.tick_bottom()

# add colourbar
ax = plt.subplot(gs[0, ncols - 1])
ax.set_xlabel("log2(Obs / Exp)\nContact Frequency")
ax.yaxis.tick_right()
plt.colorbar(img, cax=ax)
plt.savefig("Plots/apa.pileup.conditional.png")
plt.close()

# create differential heatmap
# create grid specification
ncols = 2
nrows = 1
gs = GridSpec(
    nrows=nrows,
    ncols=ncols,
    width_ratios=[20] * (ncols - 1) + [1],
)
# plotting options
plt.figure(figsize=(5 * (ncols - 1), 4 * nrows))
# make component plots
opts = dict(
    extent=[
        -(SHIFT_SIZE // 1000),  # values in kbp
        (SHIFT_SIZE // 1000),
        -(SHIFT_SIZE // 1000),
        (SHIFT_SIZE // 1000),
    ],
    cmap="bwr",
)
ax = plt.subplot(gs[0, 0])
img = ax.matshow(conditional_differential, **opts)
# add x axis labels to bottom-most subplots
ax.set_xlabel("log2(Tumour / Benign)")
ax.xaxis.set_visible(True)
ax.xaxis.tick_bottom()

# add colourbar
ax = plt.subplot(gs[:, ncols - 1])
ax.set_xlabel("log2(Obs / Exp)\nContact Frequency")
ax.yaxis.tick_right()
plt.colorbar(img, cax=ax)
plt.savefig("Plots/apa.differential.png")
plt.close()


# create differential heatmap for the top 5 and bottom 5 ranked loops
# create grid specification
ncols = len(SAMPLE_TYPES) + 3
nrows = n_top_loops
gs = GridSpec(
    nrows=nrows, ncols=ncols, width_ratios=[20] * len(SAMPLE_TYPES) + [1, 20, 1]
)
# plotting options
fig = plt.figure(figsize=(5 * (ncols - 1), 4 * nrows))
# make component plots
opts = dict(
    extent=[
        -(SHIFT_SIZE // 1000),  # values in kbp
        (SHIFT_SIZE // 1000),
        -(SHIFT_SIZE // 1000),
        (SHIFT_SIZE // 1000),
    ],
    cmap="bwr",
)
# plot benign and tumour stacks for the top and bottom loci
for i, ranked_idx in enumerate(top_loops_idx):
    for j, sample_type in enumerate(SAMPLE_TYPES):
        ax = fig.add_subplot(gs[i, j])
        img = ax.matshow(
            conditional_stack[sample_type][ranked_idx, :, :], vmin=0, vmax=4, **opts
        )
        ax.xaxis.set_visible(False)
        # add x axis labels to bottom-most subplots
        if i == nrows - 1:
            ax.set_xlabel(sample_type.title())
            ax.xaxis.set_visible(True)
            ax.xaxis.tick_bottom()

# add colourbar for conditional_stacks
ax = fig.add_subplot(gs[:, len(SAMPLE_TYPES)])
ax.set_xlabel("log2(Obs / Exp)\nContact Frequency")
ax.yaxis.tick_right()
plt.colorbar(img, cax=ax)

for i, ranked_idx in enumerate(top_loops_idx):
    # plot the difference matrix
    ax = fig.add_subplot(gs[i, len(SAMPLE_TYPES) + 1])
    img = ax.matshow(difference_stack[ranked_idx, :, :], vmin=-2, vmax=2, **opts)

# add colourbar for difference matrix
ax = fig.add_subplot(gs[:, len(SAMPLE_TYPES) + 2])
ax.set_xlabel("Tumour - Benign\nlog2(Obs / Exp)\nContact Frequency")
ax.yaxis.tick_right()
plt.colorbar(img, cax=ax)
# add x axis labels to bottom-most subplots
if i == nrows - 1:
    ax.set_xlabel("Tumour - Benign")
    ax.xaxis.set_visible(True)
    ax.xaxis.tick_bottom()


plt.savefig("Plots/apa.differential.ranked.png")
plt.close()


# create heatmap for benign-/tumour-specific loops
ncols = 3 + 1
nrows = 2
# create grid specification
gs = GridSpec(
    nrows=nrows,
    ncols=ncols,
    width_ratios=[20] * (ncols - 1) + [1],
)
# plotting options
plt.figure(figsize=(4 * (ncols - 1), 4 * nrows))
opts = dict(
    extent=[-10, 10, -10, 10],
    cmap="coolwarm",
)  # vmin=-2, vmax=2,)
# make component plots
for i, sample_type in enumerate(["tumour", "benign"]):
    for j, loop_type in enumerate(["tumour-specific", "shared", "benign-specific"]):
        ax = plt.subplot(gs[i, j])
        img = ax.matshow(np.log2(specific_apa[sample_type][loop_type]), **opts)
        # add x axis labels to bottom-most subplots
        ax.yaxis.tick_left()
        ax.xaxis.tick_bottom()
        if j == 0:
            ax.set_ylabel(sample_type)
            ax.yaxis.set_visible(True)
        if i == nrows - 1:
            ax.set_xlabel(loop_type)
            ax.xaxis.set_visible(True)


# add colourbar
ax = plt.subplot(gs[:, ncols - 1])
ax.set_xlabel("log2(Obs / Exp)\nContact Frequency")
ax.yaxis.tick_right()
plt.colorbar(img, cax=ax)
plt.savefig("Plots/apa.specific-loops.png")
plt.savefig("Plots/apa.specific-loops.pdf")
plt.close()
