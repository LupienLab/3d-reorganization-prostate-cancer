"""
test-dge
==========

Test for differential expression in samples due to connected breakpoints
"""

from __future__ import division, absolute_import, print_function
import os.path as path
import pandas as pd
import numpy as np
import networkx as nx
from interval import interval
from tqdm import tqdm
from scipy import stats
import pickle

from genomic_interval import GenomicInterval, overlapping

# ==============================================================================
# Constants
# ==============================================================================
BREAK_DIR = path.join("..", "2019-07-24_breakfinder", "Breakpoints", "Default")
TAD_DIR = path.join(
    "..", "2020-01-15_TAD-aggregation", "resolved-TADs", "separated-TADs"
)
CHRS = ["chr" + str(i) for i in list(range(1, 23)) + ["X", "Y"]]

# ==============================================================================
# Functions
# ==============================================================================


def find_tad(i: GenomicInterval, tads):
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


def get_genes_in_tads(intvls):
    """
    Query the expression table for the genes lying within the coordinates specified

    Parameters
    ----------
    intvls : Dict<str: interval.interval>
        Intervals to subset the reference genome by
    """
    # use the globally registered gene expression data
    global exprs
    # concatenate all these tables together
    involved_genes = pd.concat(
        [
            exprs.loc[
                # find genes intersecting this component for each affected component
                (exprs.chr == chrom)
                & (exprs.start <= int(component.sup))
                & (exprs.end >= component.inf),
                :,
            ]
            for chrom, intvl in intvls.items()
            for component in intvl
        ]
    )
    return involved_genes


def normalize_genes(genes, sample):
    """
    Convert expression values into a z-score for a particular sample, given the remaining

    Parameters
    ----------
    genes : pandas.DataFrame
        A set of genes to convert. Rows are genes, columns are samples
    sample : str
        ID of the sample to normalize against the remaining samples
    """
    # load the global SAMPLES variable
    global SAMPLES
    background_samples = [s for s in SAMPLES if s != sample]
    # copy this data so as to not mutate the original DataFrame
    genes_to_normalize = genes.copy()
    genes_to_normalize.set_index("EnsemblID", inplace=True)
    # calculate mean and variance of background samples for each gene
    means = genes_to_normalize[background_samples].mean(axis=1)
    # using .var instead of .std here since the background samples are the entire population, not a sample of the population
    variances = genes_to_normalize[background_samples].var(axis=1)
    # keep track of which genes have no expression in the background samples
    no_exprs_idx = variances.loc[variances == 0].index
    # keep track of which genes have no expression in the foreground sample
    sample_genes_expressed_idx = genes_to_normalize.loc[genes_to_normalize[s] > 0].index
    # calculate z-score normalization
    z = (genes_to_normalize[s] - means) / np.sqrt(variances)
    # for the samples with 0 variance, if the sample of interest also has no expression, replace NaN with 0
    z.loc[no_exprs_idx] = 0
    # for the samples with 0 variance, if the sample of interest also has some non-zero expression, replace NaN with Infinity
    z.loc[sample_genes_expressed_idx & no_exprs_idx] = np.inf
    return z


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
G = pickle.load(open("breakpoint-graphs.p", "rb"))

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


# load GENCODE reference annotation (all genes, not just protein-coding)
gencode = pd.read_csv(
    path.join("..", "..", "Data", "External", "GENCODE", "gencode.v33.all-genes.bed",),
    sep="\t",
    header=None,
    names=["chr", "start", "end", "strand", "EnsemblID"],
)

exprs = exprs.merge(gencode, on="EnsemblID", suffixes=["_exprs", "_gencode"])

# ==============================================================================
# Analysis
# ==============================================================================
# count the number of connected components in each graph
cc_counts = {s: sum(1 for cc in nx.connected_components(G[s])) for s in SAMPLES}
n_component_tests = sum(v for v in cc_counts.values())

# store hypothesis testing results
htest = pd.DataFrame.from_dict(
    {
        "Sample": np.repeat(SAMPLES, list(cc_counts.values())),
        "Component_Index": [i for n in cc_counts.values() for i in range(n)],
        "n_breakpoints": [0] * n_component_tests,
        "n_genes": [0] * n_component_tests,
        "t": [0] * n_component_tests,
        "p": [1] * n_component_tests,
    }
)
# store genes invoved in hypothesis tests
tested_genes = {
    s: {i: [] for i, _ in enumerate(nx.connected_components(G[s]))} for s in SAMPLES
}

# for each sample
# go over each component of the graph (i.e. a series of related SVs)
# and identify all the TADs that are involved (usually just 1 or 2)
# test whether the TADs in this sample are different from the TADs in the others

for s in tqdm(SAMPLES):
    for i, cc in enumerate(nx.connected_components(G[s])):
        # for each breakpoint involved in the component, find the parent TAD
        # do this across all patients and find the maximal equivalent one to get a set of sites and genes to compare
        affected_tads = pd.concat(
            [find_tad(bp, tads[new_s][3]) for new_s in SAMPLES for bp in cc]
        )
        # convert affects_tads to interval objects to more easily convert to contiguous segments later on
        affected_tad_intvls = {
            c: interval(
                *[
                    [i.start, i.end]
                    for i in affected_tads.loc[affected_tads.chr == c, :].itertuples()
                ]
            )
            for c in CHRS
        }
        # remove unused chromosomes (this is true when the value is an empty interval)
        affected_tad_intvls = {
            c: i for c, i in affected_tad_intvls.items() if i != interval()
        }
        # get all the genes across all the involved TADs
        genes = get_genes_in_tads(affected_tad_intvls)
        # CHECK THAT THIS MUTATION DOES NOT EXIST IN ANOTHER SAMPLE #
        # normalize them according to the samples without this mutation
        z = normalize_genes(genes, s)
        # order is preserved, so just append the column
        genes["z"] = z.tolist()
        # store tested genes for later
        tested_genes[s][i] = genes
        # conduct the hypothesis test, since all genes have their expression values normalized across samples
        infty_idx = genes.loc[np.isinf(genes["z"])].index
        # only conduct on non-infinite values
        # we look at these infinite values later on, specifically
        t, p = stats.ttest_1samp(genes.loc[~np.isinf(genes["z"]), "z"], 0)
        # store the data
        htest.loc[
            (htest.Sample == s) & (htest.Component_Index == i), "n_breakpoints"
        ] = len(cc)
        htest.loc[
            (htest.Sample == s) & (htest.Component_Index == i), "n_genes"
        ] = genes.shape[0] - len(infty_idx)
        htest.loc[(htest.Sample == s) & (htest.Component_Index == i), "t"] = t
        htest.loc[(htest.Sample == s) & (htest.Component_Index == i), "p"] = p

# ==============================================================================
# Save data
# ==============================================================================
htest.to_csv("sv-disruption-tests.tsv", sep="\t", index=False)
