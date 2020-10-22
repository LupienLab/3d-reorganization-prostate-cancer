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
import pandas as pd
from typing import Dict, List, Tuple
import pickle

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
metadata = pd.read_csv("config.tsv", sep="\t")
SAMPLES = {
	"all": metadata.loc[metadata.Include == "Yes", "SampleID"].tolist(),
	"tumour": metadata.loc[metadata.Type == "Malignant", "SampleID"].tolist(),
	"benign": metadata.loc[metadata.Type == "Benign", "SampleID"].tolist(),
}

# load stack data
stack = pickle.load(open(path.join(DIR["loops"], "stack.obj"), "rb"))


# ==============================================================================
# Analysis
# ==============================================================================
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

# dump to file
cdtnl_stack_obj = open(path.join(DIR["loops"], "stack.conditional.obj"), "wb")
pickle.dump(conditional_stack, cdtnl_stack_obj)
cdtnl_stack_obj.close()

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

# dump to file
difftl_stack_obj = open(path.join(DIR["loops"], "stack.conditional.differential.obj"), "wb")
pickle.dump(conditional_stack_differential, difftl_stack_obj)
difftl_stack_obj.close()

# sum over the stack of obs/exp values over each loop call
pileup = {
	s: np.nanmean(stack[s], axis=2) for s in SAMPLES["all"]
}

# dump to file
pileup_obj = open(path.join(DIR["loops"], "pileup.obj"), "wb")
pickle.dump(pileup, pileup_obj)
pileup_obj.close()

# sum over each tumour/benign sample, in groups
conditional_pileup = {
	sample_type: np.nanmean(
		[pileup[s] for s in SAMPLES[sample_type]], axis=0
	)
	for sample_type in ["tumour", "benign"]
}

# dump to file
cdntl_pileup_obj = open(path.join(DIR["loops"], "pileup.conditional.obj"), "wb")
pickle.dump(conditional_pileup, cdntl_pileup_obj)
cdntl_pileup_obj.close()

# differential testing
conditional_differential = conditional_pileup["tumour"] - conditional_pileup["benign"]

# dump to file
difftl_obj = open(path.join(DIR["loops"], "differential.obj"), "wb")
pickle.dump(conditional_differential, difftl_obj)
difftl_obj.close()

