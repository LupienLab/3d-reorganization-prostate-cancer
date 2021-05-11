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
import pandas as pd
import random
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from adjustText import adjust_text

plt.rcParams["figure.figsize"] = (20 / 2.54, 20 / 2.54)  # 20x20 cm

from genomic_interval import GenomicInterval, overlapping
import logging

import argparse

PARSER = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
PARSER.add_argument("gene_id", type=str, help="Gene ID to plot", nargs="+")
ARGS = PARSER.parse_args()


# ==============================================================================
# Constants
# ==============================================================================
# set logging parameters
logging.getLogger().setLevel(logging.INFO)

unique_annots = ["T2E-specific", "nonT2E-specific", "shared"]

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


def plot_graph(G, prefix, n_centres=10, node_labels=False, edge_labels=False, **kwargs):
    # get colours for the nodes
    n_colours = {}
    for n in G:
        if n.data["type"] == "enhancer":
            n_colours[n] = plt.cm.tab10(unique_annots.index(n.data["condition"]))
        else:
            n_colours[n] = "#000000"
    # get edges and colours corresponding to loop type
    edges, e_annots = zip(*nx.get_edge_attributes(G, "condition").items())
    e_colours = [plt.cm.tab10(unique_annots.index(a)) for a in e_annots]
    # get spring layout for the entire graph
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
        node_color=n_colours.values(),
        width=2,
        edge_color=e_colours,
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

# load gene annotations
genes = pd.read_csv(
    "../../Data/External/GENCODE/gencode.v33.all-genes.promoters.bed",
    sep="\t",
    header=None,
    names=["chr", "start", "end", "strand", "gene_id", "gene_name"],
)

# ==============================================================================
# Plots
# ==============================================================================
logging.info("Plotting GRNs")
random.seed(42)

for gid in ARGS.gene_id:
    plot_graph(
        G[gid],
        path.join(PLOT_DIR, genes.loc[genes.gene_id == gid, "gene_name"].values[0]),
        node_labels=True,
        dpi=96,
        bbox_inches="tight",
    )

# # figure for chromosome colour legend
# fig_leg, ax_leg = plt.subplots()
# cols = [locus_to_plotpos(GenomicInterval(c, 0, 1)) for c in chroms]
# ax_leg.scatter(
#     x=[a for (a, b) in cols], y=[b for (a, b) in cols], s=2000, c=chrom_colours,
# )
# ax_leg.set_ylim([0, 5])
# for i, c in enumerate(chroms):
#     ax_leg.annotate(c, cols[i], ha="center", va="center")

# plt.axis("off")
# savefig(fig_leg, path.join("Plots", "chrom-colour-map"), bbox_inches="tight")
