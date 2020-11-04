"""
plot-grn
==========

Plot GRNs for various genes
"""

import os.path as path
import networkx as nx
import pickle
from typing import List, Set, Tuple, Dict
import numpy as np
import random
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from adjustText import adjust_text

plt.rcParams["figure.figsize"] = (20 / 2.54, 20 / 2.54)  # 20x20 cm

from genomic_interval import GenomicInterval, overlapping
import logging

# ==============================================================================
# Constants
# ==============================================================================
# set logging parameters
logging.getLogger().setLevel(logging.INFO)

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

annot_map = {"INV": "V", "DEL": "Δ", "DUP": "+", "UNKNOWN": "?"}

chroms = ["chr" + str(i) for i in list(range(1, 23)) + ["X", "Y"]]

chrom_colours = [
    "#cf4336",
    "#5dbe4a",
    "#a35bce",
    "#a6b535",
    "#5e6ace",
    "#d7a036",
    "#5d88c7",
    "#d46e2b",
    "#47bcd2",
    "#dc406f",
    "#52c39a",
    "#d14da5",
    "#458631",
    "#c08ed8",
    "#7aba6f",
    "#924d89",
    "#b7ab5d",
    "#9f455c",
    "#37835c",
    "#df83a3",
    "#6a762c",
    "#bb5a53",
    "#92662d",
    "#e49670",
]

chrom_colour_map = {c: v for c, v in zip(chroms, chrom_colours)}

GRAPH_DIR = "Graphs"
PLOT_DIR = "Plots"

# ==============================================================================
# Functions
# ==============================================================================
def locus_to_plotpos(bp: GenomicInterval) -> Tuple[float, float]:
    """
    Convert genomic locus to a position in the breakpoint graphs.
    The chromosomes are laid out in a 4 rows x 6 cols pattern (chr1-22, X, Y)
    
    Parameters
    ==========
    bp: GenomicInterval
        A breakpoint to be plotted
    Returns: (float, float)
    """
    # find row and column for this breakpoint
    i = chroms.index(bp.chr)
    y = 4 - i // 6
    x = i % 6
    return (x, y)


def savefig(fig, prefix="figure", exts=["png", "pdf"], dpi=400, **kwargs):
    for ext in exts:
        fig.savefig(prefix + "." + ext, dpi=dpi, **kwargs)


def plot_graph(
    G,
    prefix,
    colocate_chroms=False,
    n_centres=10,
    node_labels=False,
    edge_labels=False,
    **kwargs
):
    if colocate_chroms:
        for c in chroms:
            for n in [m for m in G if m.chr == c]:
                for i in range(n_centres):
                    chrom_node = GenomicInterval(c, 0, 1, {"sample": None})
                    G.add_node(chrom_node)
                    G.add_edge(chrom_node, n)
        # get colours for the nodes
        chrom_nodes = [n for n in G if n.data["sample"] is None]
    # get nodes and colours corresponding to breakpoints
    bp_nodes = [n for n in G]
    # get colours for the nodes
    n_colours = [chrom_colour_map[n.chr] for n in bp_nodes]
    # get edges and colours corresponding to SV type
    edges = G.edges()
    # get spring layout for the entire graph
    if colocate_chroms:
        pos_spring = nx.spring_layout(
            G=G, pos={n: locus_to_plotpos(n) for n in chrom_nodes}, fixed=chrom_nodes,
        )
    else:
        pos_spring = nx.spring_layout(G)
    # calculate optimal distances for nodes
    optimal_dists = {}
    for n in G:
        optimal_dists[n] = {}
        for m in G:
            if G.has_edge(n, m):
                optimal_dists[n][m] = 0.2
            else:
                optimal_dists[n][m] = 1
    # use the spring_layout as the initial layout before applying the Kamada Kawai optimiziation
    pos = nx.kamada_kawai_layout(G=G, dist=optimal_dists, pos=pos_spring)
    # plot this graph
    fig, ax = plt.subplots()
    nx.draw(
        G=G,
        pos=pos,
        ax=ax,
        edgelist=edges,
        node_color=n_colours,
        linewidths=1,
        edgecolors="#000000",
        with_labels=False,
    )
    # save without node labels
    savefig(fig, prefix, **kwargs)
    if node_labels:
        labels = [0] * len(G)
        for i, n in enumerate(G):
            if n.data["type"] == "promoter":
                labels[i] = plt.text(
                    x=pos[n][0],
                    y=pos[n][1],
                    s=n.data["name"],
                    ha="center",
                    va="center",
                )
            else:
                labels[i] = plt.text(
                    x=pos[n][0], y=pos[n][1], s=n.__str__(), ha="center", va="center",
                )
        adjust_text(
            texts=labels, ax=ax, lim=10,
        )
        # save again but with node labels
        savefig(fig, prefix + ".labelled", **kwargs)
    plt.close()


# ==============================================================================
# data
# ==============================================================================
logging.info("Loading data")
# load graphs
G = pickle.load(open(path.join(GRAPH_DIR, "grns.p"), "rb"))
print(len(G))

# ==============================================================================
# Plots
# ==============================================================================
logging.info("Plotting GRNs")
random.seed(42)

counter = 0
for gene_id, grn in G.items():
    if len(grn) <= 1:
        continue
    if counter > 5:
        continue
    else:
        counter += 1
    plot_graph(
        grn,
        path.join(PLOT_DIR, gene_id),
        node_labels=True,
        dpi=96,
        bbox_inches="tight",
    )

# figure for chromosome colour legend
fig_leg, ax_leg = plt.subplots()
cols = [locus_to_plotpos(GenomicInterval(c, 0, 1)) for c in chroms]
ax_leg.scatter(
    x=[a for (a, b) in cols], y=[b for (a, b) in cols], s=2000, c=chrom_colours,
)
ax_leg.set_ylim([0, 5])
for i, c in enumerate(chroms):
    ax_leg.annotate(c, cols[i], ha="center", va="center")

plt.axis("off")
savefig(fig_leg, path.join("Plots", "chrom-colour-map"), bbox_inches="tight")
