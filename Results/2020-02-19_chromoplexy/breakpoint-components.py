"""
breakpoint-bridging
==========

Take Breakfinder results and link breakpoints together by their positions
"""

import os.path as path
from numpy import unique
import pandas as pd
import networkx as nx
from tqdm import tqdm
import pickle
from genomic_interval import GenomicInterval, overlapping

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["figure.figsize"] = (20/2.54, 20/2.54) # 20x20 cm

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
G_all = nx.Graph()
G_sample = {s: nx.Graph() for s in SAMPLES}

# add each detected breakpoint as a node to the graph
for s in tqdm(SAMPLES, unit="sample"):
    patient_bps = breakpoints.loc[s]
    for bp in patient_bps.itertuples():
        # create hashable objects to store in each node
        # in this case, and interval
        intvls = [
            GenomicInterval(
                bp.chr_x,
                bp.start_x,
                bp.end_x,
                {"sample": s, "notes": bp.notes, "strand": bp.strand_x},
            ),
            GenomicInterval(
                bp.chr_y,
                bp.start_y,
                bp.end_y,
                {"sample": s, "notes": bp.notes, "strand": bp.strand_y},
            ),
        ]
        # create nodes in the graph
        for i in intvls:
            G_all.add_node(i)
            G_sample[s].add_node(i)
        # link these two nodes since they are linked breakpoints
        G_all.add_edge(intvls[0], intvls[1], annotation=bp.annotation)
        G_sample[s].add_edge(intvls[0], intvls[1], annotation=bp.annotation)

# connect 2 nodes if their intervals are within 100 kbp of each other
for n in tqdm(G_all, position=0):
    for m in G_all:
        if n == m:
            continue
        # connect these nodes if the identified loci are within 100 kbp of each other
        if overlapping(n, m, tol / 2):
            # if the two breakpoints are from the same sample, they're nearby
            if n.data["sample"] == m.data["sample"]:
                G_all.add_edge(n, m, annotation="nearby")
            # if the breakpoints are not from the same sample, this site is recurrent
            else:
                G_all.add_edge(n, m, annotation="recurrent")

for s in SAMPLES:
    for n in tqdm(G_sample[s], position=0):
        for m in G_sample[s]:
            if n == m:
                continue
            # connect these nodes if the identified loci are within 100 kbp of each other
            if overlapping(n, m, tol / 2):
                # if the two breakpoints are from the same sample, they're nearby
                G_sample[s].add_edge(n, m, annotation="nearby")

# ==============================================================================
# Plots
# ==============================================================================
print("Plotting graphs")
# colours for edges between breakpoints
unique_annots = list(unique(breakpoints["annotation"].tolist() + ["nearby", "recurrent"]))

edges, e_annots = zip(*nx.get_edge_attributes(G_all, "annotation").items())
n_colours = [plt.cm.tab20b(SAMPLES.index(n.data["sample"])) for n in G_all]
e_colours = [plt.cm.tab10(unique_annots.index(a)) for a in e_annots]
pos = nx.circular_layout(G_all)
nx.draw(
    G_all,
    pos,
    edgelist=None,
    node_color=n_colours,
    edge_color=e_colours,
    with_labels=False,
    width=1,
)
plt.savefig(path.join("Plots", "breakpoint-graph.png"), dpi=96)
plt.close()

for i, s in enumerate(SAMPLES):
    n_colours = [plt.cm.tab20b(i)] * len(G_sample[s])
    edges, e_annots = zip(*nx.get_edge_attributes(G_sample[s], "annotation").items())
    e_colours = [plt.cm.tab10(unique_annots.index(a)) for a in e_annots]
    pos = nx.spring_layout(G_sample[s])
    nx.draw(
        G_sample[s],
        pos,
        edgelist=None,
        node_color=n_colours,
        edge_color=e_colours,
        with_labels=False,
        width=10,
    )
    plt.savefig(path.join("Plots", s + ".png"), dpi=96)
    plt.close()

# ==============================================================================
# Save data
# ==============================================================================
print("Saving graphs")
pickle.dump(G_all, open("breakpoints.all-samples.p", "wb"))
pickle.dump(G_sample, open("breakpoints.per-sample.p", "wb"))
print("Done")
