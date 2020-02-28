"""
altered-tads
==========

Determine whether SVs alter TADs
"""

from __future__ import division, absolute_import, print_function
import os.path as path
import numpy as np
import pandas as pd
import pickle
from tqdm import tqdm
from interval import interval
from genomic_interval import GenomicInterval, overlapping
from sklearn.cluster import AgglomerativeClustering

# ==============================================================================
# Constants
# ==============================================================================
BREAK_DIR = path.join("..", "2019-07-24_breakfinder", "Breakpoints", "Default")
TAD_DIR = path.join(
    "..", "2020-01-15_TAD-aggregation", "resolved-TADs", "separated-TADs"
)
CHRS = ["chr" + str(i) for i in list(range(1, 23)) + ["X", "Y"]]
WINDOWS = list(range(3, 31))
TOL = 1e5
BREAK_FLANK_SIZE = 5e5

# ==============================================================================
# Functions
# ==============================================================================
def bpscore(tads_a, tads_b, lower=None, upper=None):
    """
    Calculate the BPscore for two sets of TADs

    Parameters
    ----------
    tads_a, tads_b : 
        TADs to compare
    lower : int
        Lower bound to consider
    upper : int
        Upper bound to consider
    """
    if lower is None:
        lower = max(tads_a["start"].min(), tads_b["start"].min())
    if upper is None:
        upper = min(tads_a["end"].max(), tads_b["end"].max())
    N = upper - lower
    a = sorted(tads_a["start"].tolist() + tads_a["end"].tolist())
    # winsorize the ends so they match upper and lower
    a = np.unique([min(max(lower, x), upper) for x in a])
    b = sorted(tads_b["start"].tolist() + tads_b["end"].tolist())
    b = np.unique([min(max(lower, x), upper) for x in b])
    # counters for a, b, repsectively
    i, j = (1, 1)
    # similarity
    s = 0
    while i < len(a) and j < len(b):
        overlap = min(a[i], b[j]) - max(a[i - 1], b[j - 1])
        s += overlap ** 2 / max(a[i] - a[i-1], b[j] - b[j - 1])
        if b[j] > a[i]:
            i += 1
        else:
            j += 1
    return 1 - s / N


def different_tads(mut, mut_ids, nonmut, nonmut_ids, lower, upper):
    """
    Determine if the TADs of mutated samples are different from the TADs of
    non-mutated samples.
    Makes use of BPscore algorithm, from Zaborowski and Wilczynski 2019

    Parameters
    ----------
    mut: pd.DataFrame
        TADs of mutated samples
    mut_ids: list<str>
        List of mutated sample IDs
    nonmut: pd.DataFrame
        TADs of non-mutated samples
    nonmut_ids: list<str>
        List of non-mutated sample IDs
    lower : int
        Lower bound of positions to consider
    upper : int
        Upper bound of positions to consider
    """
    all_ids = mut_ids + nonmut_ids
    all_tads = pd.concat([mut, nonmut])
    # calculate BPscore for all pairwise comparisons
    X = pd.DataFrame(
        [[bpscore(all_tads.loc[s1], all_tads.loc[s2], lower, upper) for s2 in all_ids] for s1 in all_ids],
        index=all_ids,
        columns=all_ids
    )
    # perform kmeans clustering
    clustering = AgglomerativeClustering(n_clusters=2, affinity="precomputed", linkage="complete").fit(X)
    # compare actual mut and non-mut samples to identified clusters
    # mut IDs are always the first few rows/columns
    # thus, if they group together, the firs tlen(mut_ids) should all have the same label
    # and all the remaining IDs should have the other
    labels = clustering.labels_
    # force first label to be 1 (i.e. all mutated samples should have label 1 if the TAD changes)
    if labels[0] != 1:
        # swap labels (i.e. 0 -> 1 and 1 -> 0)
        labels = [1 - x for x in labels]
    # check consistency of clustering with mutated samples
    labels_if_TAD_altered = [1 for i in mut_ids] + [0 for i in nonmut_ids]
    if np.all(labels == labels_if_TAD_altered) :
        return True
    return False


def get_neighbouring_TADs(bp, ids, w=3, flank=BREAK_FLANK_SIZE):
    """
    Get all the TADs in a given set of samples from their IDs within a flanking distance of a given position

    Parameters
    ----------
    bp : GenomicInterval
        Breakpoint under consideration
    ids : list<str>
        Sample IDs
    w : int
        Window size parameter for TAD calls to check
    flank : int
        Flanking distance around position
    """
    # use globally loaded set of TADs
    global tads
    return pd.concat(
        [
            tads[s][w].loc[
                (tads[s][w].chr == bp.chr)
                & (tads[s][w].start <= bp.sup() + flank)
                & (bp.inf() - flank <= tads[s][w].end),
                :,
            ]
            for s in ids
        ],
        keys=ids,
    )


def get_mutated_ids_near_breakend(bp, neighbours, tol=TOL):
    """
    Get the sample IDs of all samples with a breakpoint near a given breakpoint end

    Parameters
    ----------
    bp : nx.Node
        Breakpoint under consideration
    neighbours : list<nx.Node>
        Nearby and recurrent neighbours of `bp`
    tol : int
        Distance tolerance around `pos` to be considered for "recurrent"
    """
    return set(
        [bp.data["sample"]]
        + [n.data["sample"] for n in nbrs if overlapping(bp, n, tol)]
    )


# ==============================================================================
# Data
# ==============================================================================
CONFIG = pd.read_csv(
    path.join("..", "..", "Data", "External", "LowC_Samples_Data_Available.tsv"),
    sep="\t",
    index_col=False,
)
SAMPLES = ["PCa" + str(i) for i in CONFIG["Sample ID"]]

# load TADs for each patient
tads = {
    s: {
        w: pd.read_csv(
            path.join(TAD_DIR, s + ".40000bp.w_" + str(w) + ".domains.bed"),
            sep="\t",
            names=["chr", "start", "end", "lower_persistence", "upper_persistence",],
        )
        for w in WINDOWS
    }
    for s in SAMPLES
}

# load breakpoint graphs for each patient
G = pickle.load(open("breakpoints.all-samples.p", "rb"))

# ==============================================================================
# Main
# ==============================================================================
# table to track how many breakpoint ends alter TADs
table_cols = [
    "chr",
    "start",
    "end",
    "mutated_sample",
    "altered_TAD",
]
altering_bps = pd.DataFrame(columns=table_cols)

# window size to check TADs for
w = 3

# iterate over each breakpoint
for bp in tqdm(G):
    data_to_store = {
        "chr": bp.chr,
        "start": bp.inf(),
        "end": bp.sup(),
        "mutated_sample": bp.data["sample"],
        "altered_TAD": False,
    }
    # find samples where this, or a nearby, breakpoint occurs
    nbrs = [n for n, v in G[bp].items() if v["annotation"] in ["nearby", "recurrent"]]
    # only keep unique sample IDs, don't double count them
    mut_samples = get_mutated_ids_near_breakend(bp, nbrs)
    nonmut_samples = set([s for s in SAMPLES if s not in mut_samples])
    # get TADs where this breakpoint end is in, for the mutated and non-mutated samples
    mut_tads = get_neighbouring_TADs(bp, mut_samples, w)
    nonmut_tads = get_neighbouring_TADs(bp, nonmut_samples, w)
    if different_tads(mut_tads, list(mut_samples), nonmut_tads, list(nonmut_samples), bp.inf() - BREAK_FLANK_SIZE, bp.sup() + BREAK_FLANK_SIZE):
        data_to_store["altered_TAD"] = True
    # store the result
    df = pd.DataFrame(data_to_store, index=[0])
    altering_bps = altering_bps.append(df, ignore_index=True, sort=False)

print(altering_bps["altered_TAD"].sum() / altering_bps.shape[0])

altering_bps.to_csv("sv-ends-altering-TADs.tsv", sep="\t", index=False)

