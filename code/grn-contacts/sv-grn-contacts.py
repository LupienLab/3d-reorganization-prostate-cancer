# ==============================================================================
# Meta
# ==============================================================================
# grn
# --------------------------------------
# Description: Create gene regulatory networks for each gene
# Author: James Hawley

import os.path as path
import numpy as np
import pandas as pd
import networkx as nx
from pandas.core.dtypes import dtypes
import negspy.coordinates as nc
import logging
import pickle
import itertools
from tqdm import tqdm
from typing import Tuple, Dict
import re

import sys

sys.path.insert(1, path.join("..", "generate-grns"))
from genomic_interval import GenomicInterval, overlapping, find_tad, Loop


# ==============================================================================
# Constants
# ==============================================================================
# set logging parameters
logging.getLogger().setLevel(logging.INFO)

# get chromosome sizes
hg38 = nc.get_chrominfo("hg38")
CHROM_SIZES = hg38.chrom_lengths

GRN_DIR = path.join("..", "..", "Results", "generate-grns")
SV_DIR = path.join("..", "..", "Results", "2020-02-19_chromoplexy")
EXPRS_DIR = path.join("..", "..", "Results", "2020-06-18_sv-disruption-expression")
RESULT_DIR = path.join("..", "..", "Results", "grn-contacts")


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

# load table of GRNs
G_df = pd.read_csv(
    path.join(GRN_DIR, "putative-grns.tsv"),
    sep="\t",
    header=[0],
    index_col=None,
)

# get all genes in a TAD containing an SV breakpoint for each reconstructed SV
sv_genes = pd.read_csv(
    path.join(EXPRS_DIR, "summary-sv-disruption.tsv"),
    sep="\t",
    index_col=None,
    usecols=["event_ID", "target_id", "gene_name", "status"],
)

# ==============================================================================
# Analysis
# ==============================================================================
# merge SV-related genes and their GRNs into a single table
# some genes have no putative GRNs, based on the lack of loop calls in the area
# these genes are ignored
sv_grns = sv_genes.merge(
    right=G_df,
    left_on=["target_id", "gene_name"],
    right_on=["gene_id", "gene_name"],
)

# construct table of contact frequencies between gene promoters and all enhancers involved in an SV event
eid_contacts_list = []
unique_event_IDs = list(sv_grns["event_ID"].unique())
for i, eid in tqdm(enumerate(unique_event_IDs), total=len(unique_event_IDs)):
    eid_grns = sv_grns.loc[sv_grns["event_ID"] == eid, :]
    eid_promoters = (
        eid_grns.loc[:, ["gene_id", "gene_name", "chr", "start_prom", "end_prom"]]
        .drop_duplicates()
        .reset_index(drop=True)
    )
    eid_enhancers = (
        eid_grns.loc[:, ["gene_id", "chr", "start_enh", "end_enh"]]
        .drop_duplicates()
        .reset_index(drop=True)
    )
    # match all genes with all enhancers via outer product of DataFrame indices
    prod_idx = itertools.product(eid_promoters.index, eid_enhancers.index)
    # create table of all relevant event promoter-enhancer contacts
    eid_contacts = pd.concat(
        [
            pd.DataFrame(
                {
                    "event_ID": [eid],
                    "gene_id": [eid_promoters.loc[j, "gene_id"]],
                    "gene_name": [eid_promoters.loc[j, "gene_name"]],
                    "chr_prom": [eid_promoters.loc[j, "chr"]],
                    "start_prom": [eid_promoters.loc[j, "start_prom"]],
                    "end_prom": [eid_promoters.loc[j, "end_prom"]],
                    "orig_enh_target_gene_id": [eid_enhancers.loc[k, "gene_id"]],
                    "chr_enh": [eid_enhancers.loc[k, "chr"]],
                    "start_enh": [eid_enhancers.loc[k, "start_enh"]],
                    "end_enh": [eid_enhancers.loc[k, "end_enh"]],
                }
            )
            for j, k in prod_idx
        ],
        ignore_index=True,
    )
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
contacts_all["Mutation_Status"] = contacts_all["Sample_ID"] == mut_sample
