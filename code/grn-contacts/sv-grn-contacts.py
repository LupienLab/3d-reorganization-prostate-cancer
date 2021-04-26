# ==============================================================================
# Meta
# ==============================================================================
# sv-grn-contacts
# --------------------------------------
# Description: Create gene regulatory networks for each gene
# Author: James Hawley

import os.path as path
import numpy as np
import pandas as pd
import negspy.coordinates as nc
import logging
import itertools
from tqdm import tqdm
from typing import Tuple, Dict
import pickle

from cooler import Cooler
from cooltools import snipping

import time

# ==============================================================================
# Constants
# ==============================================================================
# set logging parameters
logging.getLogger().setLevel(logging.INFO)

# get chromosome sizes
hg38 = nc.get_chrominfo("hg38")
CHROM_SIZES = hg38.chrom_lengths

# important directories
DIR = {
    "grn": path.join("..", "..", "Results", "generate-grns"),
    "sv": path.join("..", "..", "Results", "2020-02-19_chromoplexy"),
    "exprs": path.join("..", "..", "Results", "2020-06-18_sv-disruption-expression"),
    "res": path.join("..", "..", "Results", "grn-contacts"),
    "contacts": path.join(
        "..", "..", "Data", "Processed", "2019-06-18_PCa-LowC-sequencing", "Contacts"
    ),
}

# resolution for contact matrices
RES = 10000

# +/- number of bps for aggregate peak analysis
SHIFT_SIZE = 50000


# ==============================================================================
# Data
# ==============================================================================
logging.info("Loading data")

# load sample metadata
meta = pd.read_csv(
    path.join("..", "config.tsv"),
    sep="\t",
    index_col=False,
)
meta = meta.loc[(meta.Source == "Primary") & (meta.Type == "Malignant"), :]
SAMPLES = meta["Sample_ID"].tolist()

# load table of GRNs
G_df = pd.read_csv(
    path.join(DIR["grn"], "putative-grns.tsv"),
    sep="\t",
    header=[0],
    index_col=None,
)

# get all genes in a TAD containing an SV breakpoint for each reconstructed SV
sv_genes = pd.read_csv(
    path.join(DIR["exprs"], "summary-sv-disruption.tsv"),
    sep="\t",
    index_col=None,
    usecols=["event_ID", "target_id", "gene_name", "chr", "start", "end", "status"],
)
sv_genes.columns = [
    "event_ID",
    "gene_id",
    "gene_name",
    "chr",
    "start_gene",
    "end_gene",
    "gene_dge_status",
]

CHROM_SIZES = pd.read_csv(
    path.join(DIR["contacts"], "..", "hg38.sizes.txt"),
    sep="\t",
    index_col=0,
    header=None,
    names=["size"],
)
CHRS = list(CHROM_SIZES.index)

# ==============================================================================
# Analysis
# ==============================================================================
# 1. Get important loci to be extracted from contact matrices
# ------------------------------------------------
logging.info("Calculating coordinates")

# merge SV-related genes and their GRNs into a single table
# some genes have no putative GRNs, based on the lack of loop calls in the area
# these genes are ignored
sv_grns = sv_genes.merge(
    right=G_df,
    # ignore genes not relevant to SVs
    how="left",
    on=["gene_id", "gene_name", "chr", "start_gene", "end_gene"],
    validate="many_to_many",
)

# drop genes without any detected GRN (these are genes where start_enh and end_enh are NaNs)
sv_grns.dropna(axis=0, how="any", inplace=True)

# force data types for faster processing
# this didn't happen by default because of NaNs
sv_grns = sv_grns.convert_dtypes(
    {
        "start_prom": int,
        "end_prom": int,
        "start_enh": int,
        "end_enh": int,
        "detected_prom_enh_loop": bool,
    },
)
sv_grns.to_csv(
    path.join(DIR["res"], "sv-grns.tsv"),
    sep="\t",
    index=False,
    header=True,
)

# construct table of contact frequencies between gene promoters and all enhancers involved in an SV event
eid_contacts_list = []
unique_event_IDs = set(sv_grns["event_ID"])
for i, eid in tqdm(enumerate(unique_event_IDs), total=len(unique_event_IDs)):
    eid_grns = sv_grns.loc[sv_grns["event_ID"] == eid, :]
    eid_promoters = (
        eid_grns.loc[
            :, ["event_ID", "gene_id", "gene_name", "chr", "start_prom", "end_prom"]
        ]
        .drop_duplicates()
        .reset_index(drop=True)
    )
    eid_enhancers = (
        eid_grns.loc[
            :,
            ["gene_id", "chr", "start_enh", "end_enh"],
        ]
        .drop_duplicates()
        .reset_index(drop=True)
    )
    # match all genes with all enhancers via outer product of DataFrame indices
    prod_idx = list(itertools.product(eid_promoters.index, eid_enhancers.index))
    idx_proms = [idx for idx, _ in prod_idx]
    idx_enhs = [idx for _, idx in prod_idx]
    prod_proms = eid_promoters.iloc[idx_proms].reset_index(drop=True)
    prod_enhs = eid_enhancers.iloc[idx_enhs].reset_index(drop=True)
    # create table of all relevant event promoter-enhancer contacts
    eid_contacts = pd.concat(
        [prod_proms, prod_enhs],
        axis=1,
    )
    eid_contacts.columns = [
        "event_ID",
        "query_gene_id",
        "query_gene_name",
        "query_gene_prom_chr",
        "query_gene_prom_start",
        "query_gene_prom_end",
        "enh_target_gene_id",
        "enh_chr",
        "enh_start",
        "enh_end",
    ]
    eid_contacts_list.append(eid_contacts)

# combine into a single table
contacts = pd.concat(eid_contacts_list, ignore_index=True)
# expand into a table for recording the value within each locus for each sample
contacts_all = (
    pd.concat([contacts for s in meta.Sample_ID], keys=meta.Sample_ID)
    .reset_index()
    .drop("level_1", axis=1)
)

# extract mutation status for this sample from the event_IDs
mut_sample = contacts_all["event_ID"].str.extract("(PCa\d+)")[0]
contacts_all["sample_mutated_in_event"] = contacts_all["Sample_ID"] == mut_sample

# 2. Convert genomic coordinates to contact matrix indices with cooltools.snipping
# ------------------------------------------------
logging.info("Extracting contact frequencies")

# chromosomes and their sizes
supports = [(chrom, 0, CHROM_SIZES.at[chrom, "size"]) for chrom in CHRS]

# use snipping module to get chromosome locations around each gene promoter
windows_x = snipping.make_bin_aligned_windows(
    binsize=RES,
    chroms=contacts["query_gene_prom_chr"],
    # 50 kbp window centred on promoter
    centers_bp=(
        (contacts["query_gene_prom_start"] + contacts["query_gene_prom_end"]) // 2
    ),
    flank_bp=SHIFT_SIZE,
)
# do the same as above but for each enhancer
windows_y = snipping.make_bin_aligned_windows(
    binsize=RES,
    chroms=contacts["enh_chr"],
    # 50 kbp window centred on promoter
    centers_bp=(contacts["enh_start"] + contacts["enh_end"]) // 2,
    flank_bp=SHIFT_SIZE,
)

# merge the two pairs of windows together
windows_merged = pd.merge(
    left=windows_x,
    right=windows_y,
    left_index=True,
    right_index=True,
    suffixes=(
        "1",
        "2",
    ),  # these suffixes are required, assign_regions doesn't work properly if not '1' and '2'
)
# assign chromosome regions to the low/high indices
windows = snipping.assign_regions(windows_merged, supports)
windows = windows.dropna()

# 3. Extract contact frequencies from contact matrices
# ------------------------------------------------
logging.info("Calculating Obs/Exp matrices")

# contact matrices
mtx = {
    s: Cooler(path.join(DIR["contacts"], s + ".mcool::/resolutions/" + str(RES)))
    for s in SAMPLES
}

# expected matrix values
expected_mtx = {
    s: pd.read_csv(path.join(DIR["contacts"], s + ".res_10000bp.exp.cis.tsv"), sep="\t")
    for s in SAMPLES
}

# calculate obs/exp matrix for each sample and save
snipper = {s: snipping.ObsExpSnipper(mtx[s], expected_mtx[s]) for s in SAMPLES}
pickle.dump(snipper, open(path.join(DIR["res"], "snipper.obj"), "wb"))

logging.info("Aggregating matrices at loops")
# taken from https://cooltools.readthedocs.io/en/latest/notebooks/06_snipping-pileups.html
# create a stack of obs/exp matrices based on the locations in windows
# this is the part that takes the longest time
stack = {}
for s in SAMPLES:
    start_time = time.time()
    stack[s] = snipping.pileup(windows, snipper[s].select, snipper[s].snip)
    print("%s s" % round(time.time() - start_time, 2))

# save serialized object
pickle.dump(stack, open(path.join(DIR["res"], "stack.obj"), "wb"))

# calculate mean frequency over each locus
stack_mean = {s: np.nanmean(stack[s], axis=(0, 1)) for s in SAMPLES}

# merge information into large table
vals = list(itertools.chain.from_iterable([m for m in stack_mean.values()]))
contacts_all["Mean_Obs_Exp_Freq"] = vals

# ==============================================================================
# Save data
# ==============================================================================
logging.info("Saving tables")
contacts.to_csv(
    path.join(DIR["res"], "contacts.tsv"),
    sep="\t",
    index=False,
)
contacts_all.to_csv(
    path.join(DIR["res"], "contacts-extracted.tsv"),
    sep="\t",
    index=False,
)
