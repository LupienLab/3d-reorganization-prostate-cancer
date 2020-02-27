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

# ==============================================================================
# Constants
# ==============================================================================
BREAK_DIR = path.join("..", "2019-07-24_breakfinder", "Breakpoints", "Default")
TAD_DIR = path.join(
    "..", "2020-01-15_TAD-aggregation", "resolved-TADs", "separated-TADs"
)
CHRS = ["chr" + str(i) for i in list(range(1, 23)) + ["X", "Y"]]
WINDOWS = list(range(3, 31))
TOL = 100000
BREAK_FLANK_SIZE = 1000000

# ==============================================================================
# Functions
# ==============================================================================
def find_tad(i: GenomicInterval, tads):
    """
    Find the TAD(s) that overlap a given interval

    Parameters
    ----------
    i: GenomicInterval
        Interval to find TADs for
    tads: pd.DataFrame
        Set of TADs to query
    """
    # get TAD containing infimum
    tads_lower = tads.loc[
        (tads.chr == i.chr) & (tads.start <= i.inf()) & (i.inf() <= tads.end), :
    ]
    # get TAD containing supremum
    tads_upper = tads.loc[
        (tads.chr == i.chr) & (tads.start <= i.sup()) & (i.sup() <= tads.end), :
    ]
    # get indices of the TADs identified
    lower_idx = tads_lower.index.tolist()
    upper_idx = tads_upper.index.tolist()
    # if everything is contained to a single TAD, return just the one
    if (
        (len(lower_idx) == 1)
        and (len(upper_idx) == 1)
        and (lower_idx[0] == upper_idx[0])
    ):
        return tads_lower
    # if not, return the set of TADs spanned by the breakpoint coordinates
    else:
        return tads.iloc[range(min(lower_idx), max(upper_idx) + 1), :]


def different_tads(mut, nonmut):
    """
    Determine if the TADs of mutated samples are different from the TADs of
    non-mutated samples.
    Makes use of BPscore algorithm, from Zaborowski and Wilczynski 2019

    Parameters
    ----------
    mut: pd.DataFrame
        TADs of mutated samples
    nonmut: pd.DataFrame
        TADs of non-mutated samples
    """
    pass


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


def print_tads(mut, nonmut, res=40000):
    """
    Print TADs and boundaries with ASCII characters

    Parameters
    ----------
    mut : pd.DataFrame
        TADs of mutated samples
    nonmut : pd.DataFrame
        TADs of non-mutated samples
    res : int
        TAD resolution
    """
    lower = min(mut.start.min(), nonmut.start.min())
    upper = max(mut.end.max(), nonmut.end.max())
    print(mut)
    print(nonmut)
    print()
    print(lower)
    print(upper)


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
w = 10

# iterate over each breakpoint
# counter = 0
for bp in G:
    # counter += 1
    # if counter > 20:
    #     break
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
    mut_tads = get_neighbouring_TADs(bp, mut_samples)
    nonmut_tads = get_neighbouring_TADs(bp, nonmut_samples)
    print(bp)
    print(mut_tads)
    print(nonmut_tads)
    input()
    if different_tads(mut_tads, nonmut_tads):
        data_to_store[which] = True
    # store the result
    df = pd.DataFrame(data_to_store, index=[0])
    altering_bps = altering_bps.append(df, ignore_index=True, sort=False)

altering_bps["num_ends_altering_TADs"] = (
    altering_bps["lower_end_altered_TAD"] + altering_bps["upper_end_altered_TAD"]
)
print(altering_bps["num_ends_altering_TADs"].sum() / (2 * altering_bps.shape[0]))

altering_bps.to_csv("sv-ends-altering-TADs.tsv", sep="\t", index=False)

