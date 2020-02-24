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
from tqdm import tqdm
from scipy import stats

# ==============================================================================
# Constants
# ==============================================================================
BREAK_DIR = path.join("..", "2019-07-24_breakfinder", "Breakpoints", "Default")
TAD_DIR = path.join("..", "2020-01-15_TAD-aggregation", "resolved-TADs")

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
        return int(self.interval[0].inf)

    def sup(self):
        return int(self.interval[0].sup)


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


def find_tad(i: GenomicInterval, tads):
    return GenomicInterval(i.chr, tads.at[].start, tads.end)

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

# label interactions that are on different chromosomes as BND
breakpoints.loc[breakpoints.chr_x != breakpoints.chr_y, "annotation"] = "BND"

# load TADs for each patient
tads = pd.concat(
    [
        pd.read_csv(
            path.join(TAD_DIR, s + ".40000bp.aggregated-domains.sorted.bedGraph"),
            sep="\t",
            names=[
                "chr",
                "start",
                "end",
                "lower_persistence",
                "upper_persistence",
                "w",
            ],
        )
        for s in SAMPLES
    ],
    keys=SAMPLES,
)
# move sample and list index per patient to columns
tads.reset_index(inplace=True)
# rename them for easier access
tads.rename(columns={"level_0": "Sample"}, inplace=True)

# ==============================================================================
# Analysis
# ==============================================================================
# distance tolerance for comparisons
tol = 1000000

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

# store hypothesis testing results
htest = pd.DataFrame(columns=["Sample", "Component_Index", "n_breakpoints", "n_genes", "t", "p"])
# store genes invoved in hypothesis tests
tested_genes = {s: {i: [] for i, _ in enumerate(nx.connected_components(G[s]))} for s in SAMPLES}

# for each sample
# go over each component of the graph (i.e. a series of related SVs)
# and identify all the TADs that are involved (usually just 1 or 2)
# test whether the TADs in this sample are different from the TADs in the others

for s in SAMPLES:
    for i, cc in enumerate(nx.connected_components(G[s])):
        involved_tads = []
        # for each breakpoint involved in the component, find the parent TAD
        # do this across all patients and find the maximal equivalent one to get a set of sites and genes to compare
        for bp in tqdm(cc):
            tad_in_mutated_patient = find_tad(bp, tads.loc[tads.Sample == s, :])
            equivalent_tads = pd.concat([tad_in_mutated_patient] + [
                find_equivalent_tad(tad_in_mutated_patient, tads.loc[tads.Sample == new_s,:])
                for new_s in SAMPLES if new_s != s
            ])
            involved_tads.append(largest_tad(equivalent_tads))
        # get all the genes across all the involved TADs
        genes = genes_in_tad(exprs, involved_tads)
        # store tested genes for later
        tested_genes[s][i] = genes
        # CHECK THAT THIS MUTATION DOES NOT EXIST IN ANOTHER SAMPLE #
        # normalize them according to the samples without this mutation
        genes_z = normalize_genes(genes, exclude=s)
        # conduct the hypothesis test, since all genes have their expression values normalized across samples
        t, p = stats.ttest_1samp(genes_z, 0)
        # store the data
        htest.loc[(htest.Sample == s) & (htest.Component_Index == i), "n_breakpoints"] = len(cc)
        htest.loc[(htest.Sample == s) & (htest.Component_Index == i), "n_genes"] = genes_in_tad.shape[0]
        htest.loc[(htest.Sample == s) & (htest.Component_Index == i), "t"] = t
        htest.loc[(htest.Sample == s) & (htest.Component_Index == i), "p"] = p

# save data
htest.to_csv("sv-disruption-tests.tsv", sep="\t", index=False)


# ==============================================================================
# Plots 
# ==============================================================================
# colours for edges between breakpoints
unique_annots = list(np.unique(breakpoints["annotation"].tolist() + ["nearby"]))

s = "PCa57054"

edges, annots = zip(*nx.get_edge_attributes(G[s], "annotation").items())
colours = [plt.cm.tab10(unique_annots.index(a)) for a in annots]
pos = nx.spring_layout(G[s])
nx.draw(
    G[s],
    pos,
    edgelist=edges,
    edge_color=colours,
    with_labels=True,
    font_weight="bold",
    width=10,
)
