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
from genomic_interval import GenomicInterval, overlapping, find_tad

# ==============================================================================
# Constants
# ==============================================================================
BREAK_DIR = path.join("..", "2019-07-24_breakfinder", "Breakpoints", "Default")
TAD_DIR = path.join(
    "..", "2020-01-15_TAD-aggregation", "resolved-TADs", "separated-TADs"
)

# ==============================================================================
# Functions
# ==============================================================================
def equivalent_tad(
    bp_i: GenomicInterval,
    bp_j: GenomicInterval,
    tads_i: pd.DataFrame,
    tads_j: pd.DataFrame,
):
    """
    Find if two breakpoints are located in equivalent TADs in their respective samples

    Parameters
    ----------
    bp_i: GenomicInterval
        First breakpoint
    bp_j: GenomicInterval
        Second breakpoint
    tads_i: pd.DataFrame
        TADs from the sample corresponding to `bp_i`
    tads_j : pd.DataFrame
        TADs from the sample corresponding to `bp_j`
    """
    # if the breakpoints are not on the same chromosome, they can't have equivalent TADs
    if bp_i.chr != bp_j.chr:
        return False
    parent_tads_i = find_tad(bp_i, tads_i)
    parent_tads_j = find_tad(bp_j, tads_j)
    # convert to contiguous GenomicInterval to use `overlapping` function
    locus_i = GenomicInterval(
        bp_i.chr, parent_tads_i["start"].min(), parent_tads_i["end"].max()
    )
    locus_j = GenomicInterval(
        bp_j.chr, parent_tads_j["start"].min(), parent_tads_j["end"].max()
    )
    return overlapping(locus_i, locus_j)


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

# load TADs
print("Reading TADs")
tads = {
    s: {
        w: pd.read_csv(
            path.join(TAD_DIR, s + ".40000bp.w_" + str(w) + ".domains.bed"),
            sep="\t",
            names=["chr", "start", "end", "lower_persistence", "upper_persistence",],
        )
        for w in range(3, 31)
    }
    for s in SAMPLES
}


# ==============================================================================
# Analysis
# ==============================================================================
# distance tolerance for comparisons
tol = 100000

# window size to check TADs for
w = 3

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

print("Connecting proximal breakpoints")
# connect 2 nodes if their intervals are within 100 kbp of each other
for n in tqdm(G_all):
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
        # connect these nodes if the identified loci at within the equivalent TADs from their respective samples
        # (I know this isn't the most efficient way to do this, but given the number of breakpoints and samples
        # it's not that much of a concern)
        if equivalent_tad(n, m, tads[n.data["sample"]][w], tads[m.data["sample"]][w]):
            G_all.add_edge(n, m, annotation="equivalent-TAD")

for s in SAMPLES:
    for n in tqdm(G_sample[s]):
        for m in G_sample[s]:
            if n == m:
                continue
            # connect these nodes if the identified loci are within 100 kbp of each other
            if overlapping(n, m, tol / 2):
                # if the two breakpoints are from the same sample, they're nearby
                G_sample[s].add_edge(n, m, annotation="nearby")
            # connect these nodes if the identified loci at within the same TAD
            if equivalent_tad(n, m, tads[s][w], tads[s][w]):
                G_all.add_edge(n, m, annotation="same-TAD")

# ==============================================================================
# Save data
# ==============================================================================
print("Saving graphs")
pickle.dump(G_all, open("breakpoints.all-samples.p", "wb"))
pickle.dump(G_sample, open("breakpoints.per-sample.p", "wb"))
print("Done")
