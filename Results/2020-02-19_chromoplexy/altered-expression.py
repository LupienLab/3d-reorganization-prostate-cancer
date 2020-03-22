"""
altered-expression
==========

Test for differential expression in samples due to connected breakpoints
"""

import os.path as path
import pandas as pd
import numpy as np
import networkx as nx
from interval import interval
from tqdm import tqdm
from scipy import stats
import pickle

from genomic_interval import (
    GenomicInterval,
    overlapping,
    get_mutated_ids_near_breakend,
    find_tad,
)

# ==============================================================================
# Constants
# ==============================================================================
BREAK_DIR = path.join("..", "2019-07-24_breakfinder", "Breakpoints", "Default")
TAD_DIR = path.join(
    "..", "2020-01-15_TAD-aggregation", "resolved-TADs", "separated-TADs"
)
CHRS = ["chr" + str(i) for i in list(range(1, 23)) + ["X", "Y"]]

# window size to check TADs for
W = 3

# ==============================================================================
# Functions
# ==============================================================================


def genes_in_interval(intvl):
    """
    Query the expression table for the genes lying within the coordinates specified

    Parameters
    ----------
    intvl : GenomicInterval
        Interval to subset the reference genome by
    """
    # use the globally registered gene expression data
    global exprs
    return exprs.loc[
        # find genes intersecting this component for each affected component
        (exprs.chr == intvl.chr)
        & (exprs.start <= intvl.sup())
        & (exprs.end >= intvl.inf()),
        :,
    ].copy()


def normalize_genes(genes, fg_ids, bg_ids, offset=1e-3):
    """
    Convert expression values into a z-score for a particular sample, given the remaining

    Parameters
    ----------
    genes : pandas.DataFrame
        A set of genes to convert. Rows are genes, columns are samples
    fg_ids : [str]
        ID(s) of the sample(s) to bed normalized
    bg_ids : [str]
        ID(s) of the sample(s) to normalize against
    """
    # copy this data so as to not mutate the original DataFrame
    genes_to_normalize = genes.copy()
    genes_to_normalize.set_index("EnsemblID_short", inplace=True)
    # calculate mean and variance of background samples for each gene
    means = {
        "fg": genes_to_normalize.loc[:, fg_ids].mean(axis=1),
        "bg": genes_to_normalize.loc[:, bg_ids].mean(axis=1),
    }
    # sample standard deviation
    dev = genes_to_normalize[bg_ids].std(axis=1)
    # keep track of which genes have no expression in the background samples
    bg_no_exprs_idx = means["bg"].loc[means["bg"] == 0].index & dev.loc[dev == 0].index
    # keep track of which genes are expressed in the foreground sample(s)
    fg_exprs_idx = genes_to_normalize.loc[means["fg"] > 0].index
    # calculate z-score normalization
    z = (means["fg"] - means["bg"]) / dev
    # calculate fold change (log2(FPKM + offset) - log2(mean + offset) to ensure no divide by 0)
    fc = np.log2((means["fg"] + offset) / (means["bg"] + offset))
    # for the samples with 0 variance, if the sample of interest also has no expression, replace NaN with 0
    z.loc[bg_no_exprs_idx] = 0
    # for the samples with 0 variance, if the sample of interest has some non-zero expression, replace NaN with Infinity
    z.loc[fg_exprs_idx & bg_no_exprs_idx] = np.inf
    return z, fc, means, dev


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

# load TADs for each patient
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

# load breakpoint graphs for each patient
G = pickle.load(open("breakpoints.all-samples.p", "rb"))

# load gene expression data
exprs = pd.read_csv(
    path.join(
        "..",
        "..",
        "Data",
        "External",
        "CPC-GENE",
        "CPC-GENE_Chen-2019_RNAseq_rsem_gene_FPKM.13-Low-C-only.tsv",
    ),
    sep="\t",
    header=0,
)
# drop "Symbol" columns
exprs.drop("Symbol", axis=1, inplace=True)
# remove annotation version number to ensure compatibility with GENCODE
exprs["EnsemblID_short"] = exprs.EnsemblID.str.replace("\\.\\d+", "")


# load GENCODE reference annotation (all genes, not just protein-coding)
gencode = pd.read_csv(
    path.join("..", "..", "Data", "External", "GENCODE", "gencode.v33.all-genes.bed",),
    sep="\t",
    header=None,
    names=["chr", "start", "end", "strand", "EnsemblID", "name"],
)
# remove annotation version number to ensure compatibility with previous RNA-seq
gencode["EnsemblID_short"] = gencode.EnsemblID.str.replace("\\.\\d+", "")

exprs = exprs.merge(gencode, on="EnsemblID_short", suffixes=["_exprs", "_gencode"])

# re-order columns
exprs_cols = [
    "chr",
    "start",
    "end",
    "name",
    "strand",
    "EnsemblID_short",
    "EnsemblID_exprs",
    "EnsemblID_gencode",
] + SAMPLES
exprs = exprs[exprs_cols]

# ==============================================================================
# Analysis
# ==============================================================================


# number of breakpoints in the entire graph
n_bps = len(G)

# store hypothesis testing results
htest = pd.DataFrame(
    {
        "chr": [bp.chr for bp in G],
        "start": [bp.inf() for bp in G],
        "end": [bp.sup() for bp in G],
        "strand": [bp.data["strand"] for bp in G],
        "mutated_in": [""] * n_bps,
        "n_genes": [0] * n_bps,
        "n_mut": [0] * n_bps,
        "n_nonmut": [0] * n_bps,
        "t": [0.0] * n_bps,
        "p": [1.0] * n_bps,
    }
)

# store genes invoved in hypothesis tests
tested_genes = pd.DataFrame(
    columns=exprs_cols
    + [
        "z",
        "log2fold",
        "mut_mean",
        "nonmut_mean",
        "nonmut_sd",
        "mutated_in",
        "breakpoint_index",
    ]
)

# and identify all the TADs that are involved (usually just 1 or 2)
for i, bp in tqdm(enumerate(G), total=len(G)):
    # find samples where this, or a nearby, breakpoint occurs
    nbrs = [
        n
        for n, v in G[bp].items()
        if v["annotation"] in ["nearby", "recurrent", "equivalent-TAD"]
    ]
    # only keep unique sample IDs, don't double count them
    mut_samples = get_mutated_ids_near_breakend(bp, nbrs)
    nonmut_samples = [s for s in SAMPLES if s not in mut_samples]
    # find the parent TAD(s) of this breakpoint
    affected_tads = find_tad(bp, tads[bp.data["sample"]][W])
    affected_intvl = GenomicInterval(
        bp.chr, affected_tads["start"].min(), affected_tads["end"].max()
    )
    # get all the genes across all the affected interval
    genes = genes_in_interval(affected_intvl)
    # normalize them according to the samples without this mutation
    z, fc, means, dev = normalize_genes(genes, mut_samples, nonmut_samples)
    # order is preserved, so just append the column
    genes["z"] = z.tolist()
    genes["log2fold"] = fc.tolist()
    genes["mut_mean"] = means["fg"].tolist()
    genes["nonmut_mean"] = means["bg"].tolist()
    genes["nonmut_sd"] = dev.tolist()
    genes["mutated_in"] = bp.data["sample"]
    genes["n_mut"] = len(mut_samples)
    genes["n_nonmut"] = len(nonmut_samples)
    genes["breakpoint_index"] = i
    # store tested genes for later
    tested_genes = pd.concat([tested_genes, genes], ignore_index=True, sort=False)
    # conduct the hypothesis test, since all genes have their expression values normalized across samples
    finite_idx = genes.loc[~np.isinf(genes["z"])].index
    n_genes_tested = len(finite_idx)
    # only conduct on non-infinite values
    # we look at these infinite values later on, specifically
    # need at least 1 degree of freedom
    if n_genes_tested > 1:
        t, p = stats.ttest_1samp(genes.loc[finite_idx, "z"], 0)
    else:
        t = np.nan
        p = np.nan
    # store the data
    htest.loc[i, "mutated_in"] = ",".join(mut_samples)
    htest.loc[i, "n_genes"] = n_genes_tested
    htest.loc[i, "n_mut"] = len(mut_samples)
    htest.loc[i, "n_nonmut"] = len(nonmut_samples)
    htest.loc[i, "t"] = t
    htest.loc[i, "p"] = p

# ==============================================================================
# Save data
# ==============================================================================
# save hypothesis test results
htest.to_csv(
    "sv-disruption-tests.expression.tsv", sep="\t", index=True, index_label="breakpoint_index"
)

tested_genes.to_csv("sv-disruption-tests.expression.gene-level.tsv", sep="\t", index=False)
