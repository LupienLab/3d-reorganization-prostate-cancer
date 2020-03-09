"""
breakpoint-bridging
==========

Take Breakfinder results and link breakpoints together by their positions
"""

import os.path as path
import networkx as nx
import pickle
from genomic_interval import GenomicInterval, overlapping
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["figure.figsize"] = (20 / 2.54, 20 / 2.54)  # 20x20 cm

# ==============================================================================
# Constants
# ==============================================================================
BREAK_DIR = path.join("..", "2019-07-24_breakfinder", "Breakpoints", "Default")
unique_annots = [
    "nearby",
    "recurrent",
    "DUP",
    "INS",
    "INV",
    "DEL",
    "BND",
    "UNKNOWN",
    "equivalent-TAD",
]

# ==============================================================================
# Data
# ==============================================================================
# load graphs
G_all = pickle.load(open("breakpoints.all-samples.p", "rb"))
G_sample = pickle.load(open("breakpoints.per-sample.p", "rb"))
samples = list(G_sample.keys())

# ==============================================================================
# Plots
# ==============================================================================
# colours for edges between breakpoints

edges, e_annots = zip(*nx.get_edge_attributes(G_all, "annotation").items())
n_colours = [plt.cm.tab20b(samples.index(n.data["sample"])) for n in G_all]
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

for i, s in enumerate(samples):
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
        with_labels=True,
        width=10,
    )
    plt.savefig(path.join("Plots", s + ".png"), dpi=96)
    plt.close()
