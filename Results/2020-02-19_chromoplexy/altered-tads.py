"""
altered-tads
==========

Determine whether SVs alter TADs
"""

import os.path as path
import numpy as np
import pandas as pd
import pickle
from tqdm import tqdm
from interval import interval
from genomic_interval import GenomicInterval, overlapping, get_mutated_ids_near_breakend
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

# load table of tests
tests = pd.read_csv("sv-disruption-tests.tsv", sep="\t", header=0, index_col="test_index")

# load table of breakpoints
breakpoints = pd.read_csv("sv-breakpoints.tsv", sep="\t", header=0, index_col="breakpoint_index")

# ==============================================================================
# Main
# ==============================================================================
# table to track how many breakpoint ends alter TADs
table_cols = [
    "test_index",
    "breakpoint_index",
    "chr",
    "start",
    "end",
    "mutated_in",
    "altered_TAD",
]
altering_bps = pd.DataFrame(columns=table_cols)

# window size to check TADs for
w = 3

# iterate over each set of linked breakpoints
for t in tqdm(tests.itertuples(), total=tests.shape[0]):
    # extract breakpoint indices and get the breakpoints from the graph
    bp_idxs = [int(i) for i in t.breakpoint_indices.split(",")]
    bps = [n for i in bp_idxs for n in G if n.data["index"] == i]
    # get (non-)mutated sample IDs
    mut_samples = [s for s in t.mutated_in.split(",")]
    nonmut_samples = [s for s in SAMPLES if s not in mut_samples]
    # since the ends of related breakpoints can be megabases apart, we test each of them individually
    for i, bp in zip(bp_idxs, bps):
        data_to_store = {
            "test_index": t.Index,
            "breakpoint_index": i,
            "chr": bp.chr,
            "start": bp.inf(),
            "end": bp.sup(),
            "mutated_in": bp.data["sample"],
            "altered_TAD": False,
        }
        mut_tads = get_neighbouring_TADs(bp, mut_samples, w)
        nonmut_tads = get_neighbouring_TADs(bp, nonmut_samples, w)
        if different_tads(
            mut_tads,
            mut_samples,
            nonmut_tads,
            nonmut_samples,
            bp.inf() - BREAK_FLANK_SIZE,
            bp.sup() + BREAK_FLANK_SIZE
        ):
            data_to_store["altered_TAD"] = True
        altering_bps = altering_bps.append(data_to_store, ignore_index=True, sort=False)

agg_altering_bps = altering_bps.groupby("breakpoint_index").agg({"test_index": list, "chr": np.unique, "start": np.unique, "end": np.unique, "mutated_in": np.unique, "altered_TAD": list})
print(agg_altering_bps.loc[agg_altering_bps.altered_TAD.str.contains("True"), :].shape[0] / agg_altering_bps.shape[0])

agg_altering_bps.to_csv("sv-disruption-tests.TADs.tsv", sep="\t", index=True)
