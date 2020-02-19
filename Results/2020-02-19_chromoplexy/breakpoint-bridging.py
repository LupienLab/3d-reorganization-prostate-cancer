"""
breakpoint-bridging
==========

Take Breakfinder results and link breakpoints together by their positions
"""

from __future__ import division, absolute_import, print_function
import argparse
import os.path as path
import numpy as np
import pandas as pd
import networkx as nx

# ==============================================================================
# Constants
# ==============================================================================
BREAK_DIR = path.join("..", "2019-07-24_breakfinder", "Breakpoints", "Default")

# ==============================================================================
# Functions
# ==============================================================================

# ==============================================================================
# Data
# ==============================================================================
# load sample metadata
CONFIG = pd.read_csv(
    path.join("..", "..", "Data", "External", "LowC_Samples_Data_Available.tsv"),
    sep="\t",
    index_col=False,
)
SAMPLES = ["PCa" + str(i) for i in CONFIG["Sample ID"]]

# load breakpoints and concatenate tables
breakpoints = pd.concat(
    [
        pd.read_csv(
            path.join(BREAK_DIR, s + ".breaks.sorted.bedpe"),
            sep="\t",
            index_col=False,
            names=[
                "chr_x",
                "start_x",
                "end_x",
                "chr_y",
                "start_y",
                "end_y",
                "name",
                "score",
                "strand_x",
                "strand_y",
                "resolution",
            ],
        )
        for s in SAMPLES
    ],
    keys=SAMPLES,
)

# create graph
G = nx.Graph()
