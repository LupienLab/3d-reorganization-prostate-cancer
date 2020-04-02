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
import random
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from adjustText import adjust_text

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

annot_map = {
    "INV": "V",
    "DEL": "Δ",
    "DUP": "+",
    "UNKNOWN": "?"
}

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


def savefig(fig, prefix="figure", exts=["png", "pdf"], dpi=400, **kwargs):
    for ext in exts:
        fig.savefig(
            prefix + "." + ext,
            dpi=dpi,
            **kwargs
        )


def plot_graph(G, prefix, colocate_chroms=False, n_centres=10, **kwargs):
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
    bp_nodes = [n for n in G_sample[s] if n.data["sample"] is not None]
    # get colours for the nodes
    n_colours = [chrom_colour_map[n.chr] for n in bp_nodes]
    # get edges and colours corresponding to SV type
    edges, e_annots = zip(*nx.get_edge_attributes(G_sample[s], "annotation").items())
    e_colours = [plt.cm.tab10(unique_annots.index(a)) for a in e_annots]
    # get spring layout for the entire graph
    if colocate_chroms:
        pos_spring = nx.spring_layout(
            G=G,
            pos={n: locus_to_plotpos(n) for n in chrom_nodes},
            fixed=chrom_nodes,
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
    pos = nx.kamada_kawai_layout(
        G=G,
        dist=optimal_dists,
        pos=pos_spring
    )
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
    # add labels to the edges for the SV type
    nx.draw_networkx_edge_labels(
        G=G,
        pos=pos,
        ax=ax,
        edge_labels={edge: annot_map[ann] for (edge, ann) in zip(edges, e_annots) if ann in annot_map.keys()}
    )
    # save without node labels
    savefig(fig, prefix + ".no-labels", **kwargs)
    plt.close()
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
    # save again but with node labels
    savefig(fig, prefix, **kwargs)
    plt.close()


# ==============================================================================
# data
# ==============================================================================
# load graphs
G_all = pickle.load(open("breakpoints.all-samples.p", "rb"))
G_sample = pickle.load(open("breakpoints.per-sample.p", "rb"))
SAMPLES = list(G_sample.keys())

# ==============================================================================
# Plots
# ==============================================================================
plot_graph(G_all, path.join("Plots", "all-breakpoints"), dpi=400, bbox_inches="tight")

random.seed(42)
for i, s in enumerate(SAMPLES):
    print(s)
    plot_graph(G_sample[s], path.join("Plots", s), dpi=96, bbox_inches="tight")

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
savefig(fig_leg, path.join("Plots", "chrom-colour-map"), bbox_inches="tight")
