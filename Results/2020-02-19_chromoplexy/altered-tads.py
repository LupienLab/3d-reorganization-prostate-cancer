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
    non-mutated samples

    Parameters
    ----------
    mut: pd.DataFrame
        TADs of mutated samples
    nonmut: pd.DataFrame
        TADs of non-mutated samples
    """
    pass

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
    "chr", "start", "end", "mutated_sample",
    "lower_end_altered_TAD", "upper_end_altered_TAD", "lower_upper_same_TAD",
    "num_ends_altering_TADs"
] 
altering_bps = pd.DataFrame(columns=table_cols)

# iterate over each breakpoint
counter = 0
for bp in tqdm(G):
    counter += 1
    if counter > 20:
        break
    data_to_store = {
        "chr": bp.chr,
        "start": bp.inf(),
        "end": bp.sup(),
        "mutated_sample": bp.data["sample"],
        "lower_end_altered_TAD": False,
        "upper_end_altered_TAD": False,
        "lower_upper_same_TAD": False,
    }
    # find samples where this, or a nearby, breakpoint occurs
    nbrs = [n for n, v in G[bp].items() if v["annotation"] in ["nearby", "recurrent"]]
    # if they're in the same TAD, don't double count them
    if find_tad(bp, tads[bp.data["sample"]][3]).shape[0] == 1:
        data_to_store["lower_upper_same_TAD"] = True
    # for each end of the breakpoint, see if any of these neighbours has a breakpoint end within a certain distance
    for which, pos in zip(["lower_end_altered_TAD", "upper_end_altered_TAD"], [bp.inf(), bp.sup()]):
        # only keep unique sample IDs, don't double count them
        mut_samples = set([bp.data["sample"]] + [n.data["sample"] for n in nbrs if (n.chr == bp.chr) and (abs(pos - n.inf()) <= TOL or abs(pos - n.sup()) <= TOL)])
        nonmut_samples = set([s for s in SAMPLES if s not in mut_samples])
        # get TADs where this breakpoint end is in, for the mutated and non-mutated samples
        mut_tads = pd.concat(
            [
                tads[s][3].loc[
                    (tads[s][3].chr == bp.chr) & (tads[s][3].start <= pos) & (pos <= tads[s][3].end),
                    :
                ] for s in mut_samples
            ],
            keys = mut_samples
        )
        nonmut_tads = pd.concat(
            [
                tads[s][3].loc[
                    (tads[s][3].chr == bp.chr) & (tads[s][3].start <= pos) & (pos <= tads[s][3].end),
                    :
                ] for s in nonmut_samples
            ],
            keys = nonmut_samples
        )
        if different_tads(mut_tads, nonmut_tads):
            data_to_store[which] = True
    # store the result
    df = pd.DataFrame(data_to_store, index=[0])
    altering_bps = altering_bps.append(df, ignore_index=True, sort=False)

altering_bps["num_ends_altering_TADs"] = altering_bps["lower_end_altered_TAD"] + altering_bps["upper_end_altered_TAD"]
print(altering_bps["num_ends_altering_TADs"].sum() / (2*altering_bps.shape[0]))

altering_bps.to_csv("sv-ends-altering-TADs.tsv", sep="\t", index=False)


