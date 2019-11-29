"""
contact-distributions.py
==========

Calculate the distance distribution of contacts in a sample
"""

from __future__ import division, absolute_import, print_function
import argparse
import numpy as np
import pandas as pd
import cooler
from tqdm import tqdm
import seaborn as sns

# ==============================================================================
# Functions
# ==============================================================================
def count(cool):
    """
    Calculate counts of contacts at a given distance

    Parameters
    ----------
    cool : str
        Cooler file URI
    """
    c = cooler.Cooler(cool)
    # count distribution
    count_dist = pd.DataFrame(columns=["Distance", "count"], dtype=np.int32)
    # iterate through chromosomes to only collect cis contacts
    for chrom in tqdm(c.chromnames[::-1]):
        p = c.pixels().fetch(chrom)
        p["Distance"] = (p.bin2_id - p.bin1_id) * c.binsize
        # sum contact counts over distances
        p_counts = p.groupby("Distance")[["count"]].sum()
        # merge counts with previous counts
        count_dist = (
            pd.concat([count_dist, p_counts], sort=True).groupby(["Distance"]).sum()
        )
        print(count_dist)
    return count_dist


def plot_counts(count_dist, output):
    """
    Plot distribution of contact counts vs distance

    Parameters
    ----------
    count_dist : pandas.DataFrame
        1-column DataFrame with 'Distance (bp)' as the index and 'Count' and the column
    """
    g = sns.scatterplot(data=count_dist, x="Distance (bp)", y="Count")
    g.set(xscale="log", yscale="log")
    g.get_figure().savefig(output)


# ==============================================================================
# Main
# ==============================================================================


def main(cool, prefix, counts=None, verbose=0):
    """
    Main
    """
    if counts is None:
        # calculate distribution of contact counts
        if verbose > 0:
            print("Calculating contact distance distribution")
        count_dist = count(cool)
        # save as file
        count_dist.columns = ["Distance (bp)", "Count"]
        count_dist.to_csv(prefix + ".tsv", sep="\t", index=True)
    else:
        if verbose > 0:
            print("Loading count data")
        count_dist = pd.read_csv(counts, sep="\t")
    if verbose > 0:
        print("Plotting contact distance distribution")
    plot_counts(count_dist, prefix + ".png")
    if verbose > 0:
        print("Done")
    return


if __name__ == "__main__":
    PARSER = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    PARSER.add_argument("cool", type=str, help="Cooler file URI")
    PARSER.add_argument(
        "-p", "--prefix", type=str, help="Prefix for output files", default="counts"
    )
    PARSER.add_argument(
        "-c",
        "--counts",
        type=str,
        help="Skip loading counts from Cooler file, load pre-aggregated counts from this file, and only plot the data",
        default=None,
    )
    PARSER.add_argument(
        "-v", "--verbose", action="count", help="Verbosity level", default=0
    )
    ARGS = PARSER.parse_args()
    main(ARGS.cool, ARGS.prefix, ARGS.counts, ARGS.verbose)
