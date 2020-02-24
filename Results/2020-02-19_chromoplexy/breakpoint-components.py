"""
breakpoint-bridging
==========

Take Breakfinder results and link breakpoints together by their positions
"""

import os.path as path
from numpy import unique
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from tqdm import tqdm
import pickle

from genomic_interval import GenomicInterval, overlapping

# ==============================================================================
# Constants
# ==============================================================================
BREAK_DIR = path.join("..", "2019-07-24_breakfinder", "Breakpoints", "Default")


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
print("Reading breakpoints")
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

# label interactions that are on different chromosomes as BND
breakpoints.loc[breakpoints.chr_x != breakpoints.chr_y, "annotation"] = "BND"

# ==============================================================================
# Analysis
# ==============================================================================
# distance tolerance for comparisons
tol = 100000

# create graph
print("Creating graphs of breakpoints")
G = {s: nx.Graph() for s in SAMPLES}

# add each detected breakpoint as a node to the graph
for s in tqdm(SAMPLES):
    patient_bps = breakpoints.loc[s]
    for bp in patient_bps.itertuples():
        # create hashable objects to store in each node
        # in this case, and interval
        intvls = [
            GenomicInterval(
                bp.chr_x,
                bp.start_x,
                bp.end_x,
                {"notes": bp.notes, "strand": bp.strand_x},
            ),
            GenomicInterval(
                bp.chr_y,
                bp.start_y,
                bp.end_y,
                {"notes": bp.notes, "strand": bp.strand_y},
            ),
        ]
        # create nodes in the graph
        for i in intvls:
            G[s].add_node(i)
        # link these two nodes since they are linked breakpoints
        G[s].add_edge(intvls[0], intvls[1], annotation=bp.annotation)
    # connect 2 nodes if their intervals are within 100 kbp of each other
    for i, n in enumerate(G[s]):
        for j, m in enumerate(G[s]):
            if i == j:
                continue
            # connect these nodes if the identified loci are within 100 kbp of each other
            if overlapping(n, m, tol / 2):
                G[s].add_edge(n, m, annotation="nearby")

# ==============================================================================
# Plots
# ==============================================================================
print("Plotting graphs")
# colours for edges between breakpoints
unique_annots = list(unique(breakpoints["annotation"].tolist() + ["nearby"]))

for s in tqdm(SAMPLES):
    edges, annots = zip(*nx.get_edge_attributes(G[s], "annotation").items())
    colours = [plt.cm.tab10(unique_annots.index(a)) for a in annots]
    pos = nx.spring_layout(G[s])
    nx.draw(
        G[s],
        pos,
        edgelist=edges,
        edge_color=colours,
        with_labels=False,
        # font_weight="bold",
        width=10,
    )
    plt.savefig(path.join("Plots", s + ".png"))

# ==============================================================================
# Save data
# ==============================================================================
print("Saving graphs")
pickle.dump(G, open("breakpoint-graphs.p", "wb"))
print("Done")
