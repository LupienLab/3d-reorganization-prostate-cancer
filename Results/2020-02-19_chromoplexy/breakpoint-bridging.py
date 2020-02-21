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
from interval import interval
import matplotlib.pyplot as plt

# ==============================================================================
# Constants
# ==============================================================================
BREAK_DIR = path.join("..", "2019-07-24_breakfinder", "Breakpoints", "Default")

# ==============================================================================
# Classes
# ==============================================================================
class GenomicInterval:
    def __init__(self, chr, start, end, data=None):
        self.chr = chr
        self.interval = interval([start, end])
        self.data = data

    def __str__(self):
        return (
            self.chr
            + " ["
            + str(int(self.interval[0].inf))
            + ", "
            + str(int(self.interval[0].sup))
            + ")"
        )

    def __repr__(self):
        return (
            self.chr
            + " ["
            + str(int(self.interval[0].inf))
            + ", "
            + str(int(self.interval[0].sup))
            + ")"
        )

    def inf(self):
        return self.interval[0].inf

    def sup(self):
        return self.interval[0].sup


# ==============================================================================
# Functions
# ==============================================================================


def overlapping(a: GenomicInterval, b: GenomicInterval, extend: int = 0):
    """
    Boolean function to see if two GenomicIntervals overlap

    Parameters
    ----------
    a: GenomicInterval
    b: GenomicInterval
    extend: int
        Amount to extend each interval by before comparing
    """
    return (
        a.chr == b.chr
        and a.inf() - extend <= b.sup() + extend
        and b.inf() - extend <= a.sup() + extend
    )


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
            path.join(BREAK_DIR, s + ".breaks.sorted.manually-resolved.tsv"),
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
                "annotation",
                "notes",
            ],
        )
        for s in SAMPLES
    ],
    keys=SAMPLES,
)

# remove artefacts
breakpoints = breakpoints.loc[breakpoints["annotation"] != "ARTEFACT", :]

# create graph
G = {s: nx.Graph() for s in SAMPLES}

# add each detected breakpoint as a node to the graph
for s in SAMPLES:
    patient_bps = breakpoints.loc[s]
    for bp in patient_bps.itertuples():
        # create hashable objects to store in each node
        # in this case, and interval
        intvls = [
            GenomicInterval(
                bp.chr_x,
                bp.start_x,
                bp.end_x,
                {"annotation": bp.annotation, "notes": bp.notes, "strand": bp.strand_x},
            ),
            GenomicInterval(
                bp.chr_y,
                bp.start_y,
                bp.end_y,
                {"annotation": bp.annotation, "notes": bp.notes, "strand": bp.strand_y},
            ),
        ]
        # create nodes in the graph
        for i in intvls:
            G[s].add_node(i)
        # link these two nodes since they are linked breakpoints
        G[s].add_edge(intvls[0], intvls[1], colour="blue")
    # connect 2 nodes if their intervals are within 100 kbp of each other
    for i, n in enumerate(G[s]):
        for j, m in enumerate(G[s]):
            if i == j:
                continue
            # connect these nodes if the identified loci are within 100 kbp of each other
            if overlapping(n, m, 100000 / 2):
                G[s].add_edge(n, m, colour="red")


s = SAMPLES[5]
s

edges, colours = zip(*nx.get_edge_attributes(G[s], 'colour').items())
pos = nx.spring_layout(G[s])
nx.draw(G[s], pos, edgelist=edges, edge_color=colours, with_labels=True, font_weight="bold", width=10)

# these are the individual component subgraphs that are not connected to other nodes in the graph
cc = nx.connected_components(G["PCa13848"]):
