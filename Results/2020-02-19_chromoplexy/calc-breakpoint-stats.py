"""
calc-breakpoint-stats
==========

Calculate various breakpoint statistics for each sample
"""

from __future__ import division, absolute_import, print_function
import argparse
import os.path as path
import numpy as np
import pandas as pd
import networkx as nx
import pickle
from scipy import stats

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import seaborn as sns

# ==============================================================================
# Constants
# ==============================================================================
GRAPH_DIR = "Graphs"
PLOT_DIR = "Plots"

# ==============================================================================
# Functions
# ==============================================================================

# ==============================================================================
# Data
# ==============================================================================
# load sample metadata
metadata = pd.read_csv(
    path.join("..", "..", "Data", "External", "LowC_Samples_Data_Available.tsv"),
    sep="\t",
    header=0
)
metadata["Sample ID"] = ["PCa" + str(i) for i in metadata["Sample ID"]]
SAMPLES = metadata["Sample ID"].tolist()

# load breakpoint graphs for each sample
G_sample = pickle.load(open(
    path.join(GRAPH_DIR, "breakpoints.per-sample.merged-breakpoints.p"),
    "rb"
))

# ==============================================================================
# Analysis
# ==============================================================================
# calculate the number of complex events in each sample
complex_events = metadata.loc[:, ["Sample ID", "T2E Status"]]
complex_events["Simple Events"] = [
    sum(1 for cc in nx.connected_components(G_sample[s]) if len(cc) == 2) for s in SAMPLES
]
complex_events["Complex Events"] = [
    sum(1 for cc in nx.connected_components(G_sample[s]) if len(cc) > 2) for s in SAMPLES
]

# calculate Mann-Whitney U test for differences in number of complex events between T2E and non-T2E samples
mwu = stats.mannwhitneyu(
    x = complex_events.loc[complex_events["T2E Status"] == "Yes", "Complex Events"],
    y = complex_events.loc[complex_events["T2E Status"] == "No", "Complex Events"],
    alternative="greater"
)

all_events = pd.concat([
    pd.DataFrame({
        "Sample ID": [s] * sum(1 for cc in nx.connected_components(G_sample[s])),
        "Length": [len(cc) for cc in nx.connected_components(G_sample[s])]
    }) for s in SAMPLES
])

event_counts = sns.catplot(
    x="Length",
    kind="count",
    hue="Sample ID",
    edgecolor="#000000",
    data=all_events
)
event_counts.savefig(
    path.join(PLOT_DIR, "event-counts.png"),
    dpi=400,
    bbox_inches="tight",
)
