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
from adjustText import adjust_text
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
    y = 4 - i // 6
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
for i, s in enumerate(samples):
    print(s)
    # get nodes and colours corresponding to breakpoints
    bp_nodes = [n for n in G_sample[s] if n.data["sample"] is not None]
    # get colours for the nodes
    #n_colours = [plt.cm.tab20b(i)] * len(bp_nodes)
    n_colours = [chrom_colour_map[n.chr] for n in bp_nodes]
    # get edges and colours corresponding to SV type
    edges, e_annots = zip(*nx.get_edge_attributes(G_sample[s], "annotation").items())
    e_colours = [plt.cm.tab10(unique_annots.index(a)) for a in e_annots]
    # get spring layout for the entire graph
    pos = nx.spring_layout(G_sample[s])
    # plot this graph
    fig, ax = plt.subplots()
    nx.draw(
        G=G_sample[s],
        pos=pos,
        ax=ax,
        edgelist=edges,
        node_color=n_colours,
        #edge_color=e_colours,
        with_labels=False,
        #width=10,
    )
    fig.tight_layout()
    fig.savefig(path.join("Plots", s + ".no-labels.png"), dpi=96, bbox_inches="tight")
    fig.savefig(path.join("Plots", s + ".no-labels.pdf"), dpi=96, bbox_inches="tight")
    labels = [
        plt.text(
            x=pos[n][0],
            y=pos[n][1],
            s=n.__str__(),
            ha="center",
            va="center",
        ) for n in G_sample[s]
    ]
    adjust_text(
        texts=labels,
        ax=ax,
        lim=10,
    )
    fig.savefig(path.join("Plots", s + ".png"), dpi=96, bbox_inches="tight")
    fig.savefig(path.join("Plots", s + ".pdf"), dpi=96, bbox_inches="tight")
    plt.close()

# figure for chromosome colour legend
fig_leg, ax_leg = plt.subplots()
cols = [locus_to_plotpos(GenomicInterval(c, 0, 1)) for c in chroms]
ax_leg.scatter(
    x=[a for (a, b) in cols],
    y=[b for (a, b) in cols],
    s=2000,
    c=chrom_colours,
)
ax_leg.set_ylim([0, 5])
for i, c in enumerate(chroms):
    ax_leg.annotate(c, cols[i], ha="center", va="center")

plt.axis("off")
fig_leg.savefig(path.join("Plots", "chrom-colour-map.png"))
