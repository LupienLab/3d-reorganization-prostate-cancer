"""
breakpoint-bridging
==========

Take Breakfinder results and link breakpoints together by their positions
"""

import os.path as path
import networkx as nx
import pickle
import numpy as np
from genomic_interval import GenomicInterval, overlapping
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["figure.figsize"] = (20 / 2.54, 20 / 2.54)  # 20x20 cm

# ==============================================================================
# Constants
# ==============================================================================
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
chroms = ["chr" + str(i) for i in list(range(1, 23)) + ["X", "Y"]]

# ==============================================================================
# Functions
# ==============================================================================
def locus_to_plotpos(bp: GenomicInterval) -> (float, float):
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
    y = 6 - i // 6
    x = i % 6
    return (x, y)

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
# colours for edges between breakpoints (only get edges with an annotation, not the ones made for positioning)
#edges, e_annots = zip(*nx.get_edge_attributes(G_all, "annotation").items())
#n_colours = [plt.cm.tab20b(samples.index(n.data["sample"])) for n in G_all]
#e_colours = [plt.cm.tab10(unique_annots.index(a)) for a in e_annots]
#pos = nx.circular_layout(G_all)
#nx.draw(
#    G_all,
#    pos,
#    edgelist=None,
#    node_color=n_colours,
#    edge_color=e_colours,
#    with_labels=False,
#    width=1,
#)
#plt.savefig(path.join("Plots", "breakpoint-graph.png"), dpi=96)
#plt.close()

# create phantom points for each chromosome and place them in a grid
# this, coupled with the edges and spring_layout will approximately force
# breakpoints on a particular chromosome to roughly be in the same spot
# these points and edges are then not shown in the resulting graph
# create weights between nodes based on their positions just for plotting
n_centres = 4
for i, s in enumerate(samples):
    # create phantom points for each chromosome and place them in a grid
    # this, coupled with the edges and spring_layout will approximately force
    # breakpoints on a particular chromosome to roughly be in the same spot
    # these points and edges are then not shown in the resulting graph
    for c in chroms:
        for j in range(n_centres):
            chrom_node = GenomicInterval(c, 0, 1, {"sample": None})
            G_sample[s].add_node(chrom_node)

for i, s in enumerate(samples):
    # get nodes and colours corresponding to breakpoints
    bp_nodes = [n for n in G_sample[s] if n.data["sample"] is not None]
    # get colours for the nodes
    n_colours = [plt.cm.tab20b(i)] * len(bp_nodes)
    # get edges and colours corresponding to SV type
    edges, e_annots = zip(*nx.get_edge_attributes(G_sample[s], "annotation").items())
    e_colours = [plt.cm.tab10(unique_annots.index(a)) for a in e_annots]
    # get phantom nodes for each chromosomes
    chrom_nodes = [n for n in G_sample[s] if n.data["sample"] is None]
    # add edges between each breakpoint and its chromosome's phantom node
    for n in bp_nodes:
        for cn in [m for m in chrom_nodes if m.chr == n.chr]:
            G_sample[s].add_edge(n, cn)
    # get spring layout for the entire graph, fixing the phantom nodes at pre-defined points
    pos = nx.spring_layout(
        G=G_sample[s],
        fixed=chrom_nodes,
        pos={n: locus_to_plotpos(n) for n in chrom_nodes},
    )
    # get positions for bp_nodes, not the phantom ones
    nodes_pos = {k: v for k, v in pos.items() if k in bp_nodes}
    # plot this graph
    fig, ax = plt.subplots()
    nx.draw(
        G_sample[s].subgraph(bp_nodes),
        nodes_pos,
        edgelist=edges,
        node_color=n_colours,
        #edge_color=e_colours,
        with_labels=True,
        #width=10,
    )
    #fig.tight_layout()
    fig.savefig(path.join("Plots", s + ".png"), dpi=96, bbox_inches="tight")
    plt.close()
