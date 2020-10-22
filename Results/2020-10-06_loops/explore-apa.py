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
from scipy.stats import ttest_1samp

# ==============================================================================
# Constants
# ==============================================================================
logging.getLogger().setLevel(logging.INFO)

DIR = {"loops": "Loops"}

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
stack = pickle.load(open(path.join(DIR["loops"], "stack.obj"), "rb"))

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
conditional_stack_differential = np.array([
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
])

# get index of loops sorted by their desired ranking (most desired at the beginning of the list)
conditional_stack_ranking = np.argsort(-conditional_stack_differential)


# ==============================================================================
# Analysis
# ==============================================================================
# sum over the stack of obs/exp values over each loop call
pileup = {
	s: np.nanmean(stack[s], axis=2) for s in SAMPLES["all"]
}

# sum over each tumour/benign sample, in groups
conditional_pileup = {
	sample_type: np.nanmean(
		[pileup[s] for s in SAMPLES[sample_type]], axis=0
	)
	for sample_type in ["tumour", "benign"]
}

# differential testing
conditional_differential = np.log2(
	conditional_pileup["tumour"]
	/ conditional_pileup["benign"]
)
differential_test = ttest_1samp(
	conditional_differential.flatten(),  # flatten into array instead of matrix
	popmean=0,
)

# ==============================================================================
# Plots
# ==============================================================================
# create heatmap for each stack plot
# get rows and columns per loop type
ncols = len(SAMPLES["all"]) + 1
nrows = np.min([stack["PCa13266"].shape[2], 20])
# create grid specification
gs = GridSpec(
	nrows=nrows, ncols=ncols, width_ratios=[20] * len(SAMPLES["all"]) + [1],
)
# plotting options
plt.figure(figsize=(4 * (ncols - 1), 4 * nrows))
opts = dict(extent=[-10, 10, -10, 10], cmap="coolwarm", vmin=-2, vmax=2,)
# make component plots
for i, ranked_idx in enumerate(conditional_stack_ranking[0:nrows]):
	for j, s in enumerate(SAMPLES["all"]):
		ax = plt.subplot(gs[i, j])
		img = ax.matshow(np.log2(stack[s][:, :, ranked_idx]), **opts)
		ax.xaxis.set_visible(False)
		if j > 0:
			ax.yaxis.set_visible(False)
# add colourbar
ax = plt.subplot(gs[:, ncols - 1])
ax.set_xlabel("log2(Obs / Exp)\nContact Frequency")
ax.yaxis.tick_right()
plt.colorbar(img, cax=ax)
plt.savefig("Plots/apa.stack.loops-by-samples.png")
plt.close()


# create heatmap for each pileup plot
# get rows and columns per loop type
ncols = len(SAMPLES["all"]) + 1
nrows = 1
# create grid specification
gs = GridSpec(nrows=nrows, ncols=ncols, width_ratios=[20] * (ncols - 1) + [1],)
# plotting options
plt.figure(figsize=(4 * (ncols - 1), 4 * nrows))
opts = dict(extent=[-10, 10, -10, 10], cmap="coolwarm",)  # vmin=-2, vmax=2,)
# make component plots
for j, s in enumerate(SAMPLES["all"]):
	ax = plt.subplot(gs[0, j])
	img = ax.matshow(np.log2(pileup[s]), **opts)
	ax.xaxis.set_visible(False)
	# add x axis labels to bottom-most subplots
	if i == nrows - 1:
		ax.set_xlabel(s)
		ax.xaxis.set_visible(True)
		ax.xaxis.tick_bottom()


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
nrows = 1
# create grid specification
gs = GridSpec(nrows=nrows, ncols=ncols, width_ratios=[20] * (ncols - 1) + [1],)
# plotting options
plt.figure(figsize=(4 * (ncols - 1), 4 * nrows))
colourvals = {"min": 0.6, "max": 1.8,}
# make component plots
opts = dict(
	extent=[-10, 10, -10, 10],
	cmap="coolwarm",
	vmin=colourvals["min"],
	vmax=colourvals["max"],
)
for j, sample_type in enumerate(SAMPLE_TYPES):
	ax = plt.subplot(gs[0, j])
	img = ax.matshow(np.log2(conditional_pileup[sample_type]), **opts)
	ax.xaxis.set_visible(False)
	# add x axis labels to bottom-most subplots
	if i == nrows - 1:
		ax.set_xlabel(sample_type.title())
		ax.xaxis.set_visible(True)
		ax.xaxis.tick_bottom()
# add colourbar
ax = plt.subplot(gs[0, ncols - 1])
ax.set_xlabel("log2(Obs / Exp)\nContact Frequency")
ax.yaxis.tick_right()
plt.colorbar(img, cax=ax)


plt.savefig("Plots/apa.condition-pileup.loop-type-by-sample-type.png")
plt.close()


# create differential heatmap for each conditional pileup plot
# get rows and columns per loop type
ncols = 2
nrows = 1
# create grid specification
gs = GridSpec(nrows=nrows, ncols=ncols, width_ratios=[20] * (ncols - 1) + [1],)
# plotting options
plt.figure(figsize=(5 * (ncols - 1), 4 * nrows))
colourvals = {"min": -0.6, "max": 0.6,}
# make component plots
opts = dict(
	extent=[-25, 25, -25, 25],
	cmap="bwr",
	vmin=colourvals["min"],
	vmax=colourvals["max"],
)
ax = plt.subplot(gs[0, 0])
img = ax.matshow(conditional_differential, **opts)
ax.xaxis.set_visible(False)
# add x axis labels to bottom-most subplots
if i == nrows - 1:
	ax.set_xlabel("log2(Tumour / Benign)")
	ax.xaxis.set_visible(True)
	ax.xaxis.tick_bottom()

# add colourbar
ax = plt.subplot(gs[:, ncols - 1])
ax.set_xlabel("log2(Obs / Exp)\nContact Frequency")
ax.yaxis.tick_right()
plt.colorbar(img, cax=ax)


plt.savefig("Plots/apa.condition-pileup.loop-type-differential.png")
plt.close()
