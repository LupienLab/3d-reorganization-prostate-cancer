"""
broken-loops
==========

Measure how the loops are intercepted by SV breakpoints
"""

from __future__ import division, absolute_import, print_function
import os.path as path
import numpy as np
import pandas as pd
import networkx as nx
import pickle


# ==============================================================================
# Main
# ==============================================================================
# load sample metadata
metadata = pd.read_csv(
    path.join("..", "..", "Data", "External", "LowC_Samples_Data_Available.tsv"),
    sep="\t",
    header=[0],
)
SAMPLES = ["PCa" + str(i) for i in CONFIG.loc[CONFIG.Include == "Yes", "Sample ID"]]

# load SVs
breakpoints = pickle.load(open("breakpoints.per-sample.p", "rb"))


# load loops
loops = pd.concat(
    [
        pd.read_csv(
            path.join(
                "..",
                "2020-01-02_loops",
                "Classification",
                s + ".loops.tsv",
            ),
            sep="\t",
            index_col=False,
        )
        for s in SAMPLES
    ],
    keys=SAMPLES
)
