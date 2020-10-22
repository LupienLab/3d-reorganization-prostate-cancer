# ==============================================================================
# Meta
# ==============================================================================
# Explore APA
# --------------------------------------
# Description: Explore APA results
# Author: James Hawley

import os.path as path
import numpy as np
import pandas as pd
import pickle

import matplotlib as mpl

mpl.use("agg")

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import logging
from datatable import dt, f, by, update
from typing import Dict, List, Tuple
from tqdm import tqdm
from scipy.stats import ttest_1samp

# ==============================================================================
# Constants
# ==============================================================================
logging.getLogger().setLevel(logging.INFO)

DIR = {"loops": "Loops"}

LOOP_TYPES = ["shared", "benign", "tumour"]
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
metadata = dt.fread("config.tsv", sep="\t")
SAMPLES = {
    "all": metadata[f.Include == "Yes", f.SampleID].to_list()[0],
    "tumour": metadata[f.Type == "Malignant", "SampleID"].to_list()[0],
    "benign": metadata[f.Type == "Benign", "SampleID"].to_list()[0],
}

# load stack data
stack = {
    loop_type: pickle.load(
        open(path.join(DIR["loops"], ".".join(["stack", loop_type, "obj"])), "rb")
    )
    for loop_type in LOOP_TYPES
}  # type: Dict[str, Dict[str, np.ndarray]]

# sum over groups for each locus to rank the differential loop enrichment across loci
conditional_stack = {
    loop_type: {
        # 3. coerce into a 3D array: (loop locus, row, col)
        sample_type: np.array(
            [
                # 2. get the sample-type mean of the obs/exp matrices for the given locus
                np.nanmean(
                    # 1. get all obs/exp matrices from samples of the specific type for a given locus, i
                    [stack[loop_type][s][:, :, i] for s in SAMPLES[sample_type]],
                    axis=0,
                )
                for i in range(stack[loop_type]["PCa13266"].shape[2])
            ]
        )
        for sample_type in SAMPLE_TYPES
    }
    for loop_type in LOOP_TYPES
}  # type: Dict[str, Dict[str, np.ndarray]]

# get differential obs/exp matrices for each loop locus
conditional_stack_differential = {
    loop_type:
    # 3. coerce into a 3D array: (loop locus, row, col)
    np.array(
        [
            # 1. calculate differential between tumour and benign for a given locus, i
            np.nanmean(
                filter_arr(
                    np.log2(
                        conditional_stack[loop_type]["tumour"][i, :, :]
                        / conditional_stack[loop_type]["benign"][i, :, :]
                    )
                )
            )
            for i in range(stack[loop_type]["PCa13266"].shape[2])
        ]
    )
    for loop_type in LOOP_TYPES
}  # type: Dict[str, Dict[str, np.ndarray]]

# get index of loops sorted by their desired ranking (most desired at the beginning of the list)
conditional_stack_ranking = {
    # smallest differences between tumour and benign
    "shared": np.argsort(1 / np.absolute(conditional_stack_differential["shared"])),
    # greatest enrichment in benign
    "benign": np.argsort(conditional_stack_differential["benign"]),
    "tumour": np.argsort(-conditional_stack_differential["tumour"]),
}


# ==============================================================================
# Analysis
# ==============================================================================
# sum over the stack of obs/exp values over each loop call
pileup = {
    loop_type: {s: np.nanmean(stack[loop_type][s], axis=2) for s in SAMPLES["all"]}
    for loop_type in LOOP_TYPES
}  # type: Dict[str, Dict[str, np.ndarray]]

# sum over each tumour/benign sample, in groups
conditional_pileup = {
    (loop_type, sample_type): np.nanmean(
        [pileup[loop_type][s] for s in SAMPLES[sample_type]], axis=0
    )
    for loop_type in LOOP_TYPES
    for sample_type in ["tumour", "benign"]
}  # type: Dict[Tuple[str, str], np.ndarray]

# differential testing
conditional_differential = {
    loop_type: np.log2(
        conditional_pileup[(loop_type, "tumour")]
        / conditional_pileup[(loop_type, "benign")]
    )
    for loop_type in LOOP_TYPES
}

differential_tests = {
    loop_type: ttest_1samp(
        conditional_differential[
            loop_type
        ].flatten(),  # flatten into array instead of matrix
        popmean=0,
    )
    for loop_type in LOOP_TYPES
}

# ==============================================================================
# Plots
# ==============================================================================
# create heatmap for each stack plot
for loop_type in LOOP_TYPES:
    # get rows and columns per loop type
    ncols = len(SAMPLES["all"]) + 1
    nrows = np.min([stack[loop_type]["PCa13266"].shape[2], 20])
    # create grid specification
    gs = GridSpec(
        nrows=nrows, ncols=ncols, width_ratios=[20] * len(SAMPLES["all"]) + [1],
    )
    # plotting options
    plt.figure(figsize=(4 * (ncols - 1), 4 * nrows))
    opts = dict(extent=[-10, 10, -10, 10], cmap="coolwarm", vmin=-2, vmax=2,)
    # make component plots
    for i, ranked_idx in enumerate(conditional_stack_ranking[loop_type][0:nrows]):
        for j, s in enumerate(SAMPLES["all"]):
            ax = plt.subplot(gs[i, j])
            img = ax.matshow(np.log2(stack[loop_type][s][:, :, ranked_idx]), **opts)
            ax.xaxis.set_visible(False)
            if j > 0:
                ax.yaxis.set_visible(False)
    # add colourbar
    ax = plt.subplot(gs[:, ncols - 1])
    ax.set_xlabel("log2(Obs / Exp)\nContact Frequency")
    ax.yaxis.tick_right()
    plt.colorbar(img, cax=ax)
    plt.savefig("Plots/apa.stack.loops-by-samples.{}.png".format(loop_type))
    plt.close()


# create heatmap for each pileup plot
# get rows and columns per loop type
ncols = len(SAMPLES["all"]) + 1
nrows = len(LOOP_TYPES)
# create grid specification
gs = GridSpec(nrows=nrows, ncols=ncols, width_ratios=[20] * (ncols - 1) + [1],)
# plotting options
plt.figure(figsize=(4 * (ncols - 1), 4 * nrows))
opts = dict(extent=[-10, 10, -10, 10], cmap="coolwarm",)  # vmin=-2, vmax=2,)
# make component plots
for i, loop_type in enumerate(LOOP_TYPES):
    for j, s in enumerate(SAMPLES["all"]):
        ax = plt.subplot(gs[i, j])
        img = ax.matshow(np.log2(pileup[loop_type][s]), **opts)
        ax.xaxis.set_visible(False)
        # add x axis labels to bottom-most subplots
        if i == nrows - 1:
            ax.set_xlabel(s)
            ax.xaxis.set_visible(True)
            ax.xaxis.tick_bottom()
        # add y axis labels to left-most subplots
        if j == 0:
            ax.set_ylabel(loop_type.title())


# add colourbar
ax = plt.subplot(gs[:, ncols - 1])
ax.set_xlabel("log2(Obs / Exp)\nContact Frequency")
ax.yaxis.tick_right()
plt.colorbar(img, cax=ax)
plt.savefig("Plots/apa.pileup.loop-type-by-samples.png")
plt.close()


# create heatmap for each conditional pileup plot
# get rows and columns per loop type
ncols = len(SAMPLE_TYPES) + 1
nrows = len(LOOP_TYPES)
# create grid specification
gs = GridSpec(nrows=nrows, ncols=ncols, width_ratios=[20] * (ncols - 1) + [1],)
# plotting options
plt.figure(figsize=(4 * (ncols - 1), 4 * nrows))
colourvals = {
    "shared": {"min": 0.6, "max": 1.8,},
    "benign": {"min": 0.3, "max": 1.4,},
    "tumour": {"min": 0.6, "max": 1.5,},
}
# make component plots
for i, loop_type in enumerate(LOOP_TYPES):
    opts = dict(
        extent=[-10, 10, -10, 10],
        cmap="coolwarm",
        vmin=colourvals[loop_type]["min"],
        vmax=colourvals[loop_type]["max"],
    )
    for j, sample_type in enumerate(SAMPLE_TYPES):
        ax = plt.subplot(gs[i, j])
        img = ax.matshow(np.log2(conditional_pileup[(loop_type, sample_type)]), **opts)
        ax.xaxis.set_visible(False)
        # add x axis labels to bottom-most subplots
        if i == nrows - 1:
            ax.set_xlabel(sample_type.title())
            ax.xaxis.set_visible(True)
            ax.xaxis.tick_bottom()
        # add y axis labels to left-most subplots
        if j == 0:
            ax.set_ylabel(loop_type.title())
    # add colourbar
    ax = plt.subplot(gs[i, ncols - 1])
    ax.set_xlabel("log2(Obs / Exp)\nContact Frequency")
    ax.yaxis.tick_right()
    plt.colorbar(img, cax=ax)


plt.savefig("Plots/apa.condition-pileup.loop-type-by-sample-type.png")
plt.close()


# create differential heatmap for each conditional pileup plot
# get rows and columns per loop type
ncols = 2
nrows = len(LOOP_TYPES)
# create grid specification
gs = GridSpec(nrows=nrows, ncols=ncols, width_ratios=[20] * (ncols - 1) + [1],)
# plotting options
plt.figure(figsize=(5 * (ncols - 1), 4 * nrows))
colourvals = {
    "shared": {"min": -0.6, "max": 0.6,},
    "benign": {"min": -0.6, "max": 0.6,},
    "tumour": {"min": -0.6, "max": 0.6,},
}
# make component plots
for i, loop_type in enumerate(LOOP_TYPES):
    opts = dict(
        extent=[-25, 25, -25, 25],
        cmap="bwr",
        vmin=colourvals[loop_type]["min"],
        vmax=colourvals[loop_type]["max"],
    )
    ax = plt.subplot(gs[i, 0])
    img = ax.matshow(conditional_differential[loop_type], **opts)
    ax.xaxis.set_visible(False)
    # add x axis labels to bottom-most subplots
    if i == nrows - 1:
        ax.set_xlabel("log2(Tumour / Benign)")
        ax.xaxis.set_visible(True)
        ax.xaxis.tick_bottom()
    # add y axis labels to left-most subplots
    ax.set_ylabel(loop_type.title())

# add colourbar
ax = plt.subplot(gs[:, ncols - 1])
ax.set_xlabel("log2(Obs / Exp)\nContact Frequency")
ax.yaxis.tick_right()
plt.colorbar(img, cax=ax)


plt.savefig("Plots/apa.condition-pileup.loop-type-differential.png")
plt.close()
