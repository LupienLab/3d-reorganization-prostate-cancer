"""
bpscore
==========

Calculate the BPscore for a set of samples
"""

import numpy as np
import pandas as pd

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
        s += overlap ** 2 / max(a[i] - a[i - 1], b[j] - b[j - 1])
        if b[j] > a[i]:
            i += 1
        else:
            j += 1
    return 1 - s / N

