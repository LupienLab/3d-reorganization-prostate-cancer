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
from tqdm import tqdm
from typing import Tuple, Dict

from genomic_interval import GenomicInterval, overlapping, find_tad, Loop

# ==============================================================================
# Constants
# ==============================================================================
# set logging parameters
logging.getLogger().setLevel(logging.INFO)

# offset distance for detecting overlaps with loop calls
LOOP_OFFSET = 5000

# offset distance of detecting loops in the first place (Mustache doesn't call loops < 20 kbp, by default)
LOOP_DISTANCE_THRESH = 20000

# get chromosome sizes
hg38 = nc.get_chrominfo("hg38")
CHROM_SIZES = hg38.chrom_lengths

RESULT_DIR = path.join("..", "..", "Results", "generate-grns")
TAD_DIR = path.join(
    "..",
    "..",
    "Results",
    "2020-08-29_TADs-downsampled",
    "Aggregated-TADs",
    "separated-TADs",
)


# ==============================================================================
# Data
# ==============================================================================
logging.info("Loading data")

# load pre-processed genes
genes = pd.read_csv(
    path.join(RESULT_DIR, "overlapping.genes.tsv"),
    sep="\t",
    header=[0],
)

# load catalogue of H3K27ac peaks and differential test results
enhancers = pd.read_csv(
    path.join(RESULT_DIR, "overlapping.enhancers.tsv"),
    sep="\t",
    header=[0],
)

# load catalogue of loop calls from all 17 samples
loops = pd.read_csv(
    path.join(RESULT_DIR, "overlapping.loops.tsv"),
    sep="\t",
    header=[0],
)

# load sample metadata
meta = pd.read_csv(
    "config.tsv",
    sep="\t",
    header=[0],
)
SAMPLES = meta.loc[meta.Include == "Yes", "Sample_ID"].tolist()

# load TAD calls
tads = pd.concat(
    [
        pd.read_csv(
            path.join(
                TAD_DIR,
                s + ".300000000.res_40000bp.window_20.domains.tsv",
            ),
            sep="\t",
            header=None,
            names=[
                "chr",
                "start",
                "end",
                "persistence_left",
                "persistence_right",
                "type",
            ],
        )
        for s in SAMPLES
    ],
    keys=SAMPLES,
)
tads.reset_index(inplace=True)
tads.rename({"level_0": "Sample_ID"}, axis=1, inplace=True)
tads.drop(labels=["level_1"], axis=1, inplace=True)
# sort by chr, start, end, SampleID
tads.sort_values(
    by=["chr", "start", "end", "Sample_ID"], inplace=True, ignore_index=True
)


# ==============================================================================
# Analysis
# ==============================================================================
logging.info("Creating initial GRNs")

# create placeholder for the GRN
G = {gid: nx.Graph() for gid in genes["gene_id"]}

T = {}  # TADs
P = {}  # Genes
for gene in tqdm(genes.itertuples(), total=genes.shape[0]):
    # restrict to just the promoter for the GRN
    gi = GenomicInterval(
        gene.chr,
        gene.start - 1500 if gene.strand == "+" else gene.end + 1500,
        gene.start + 500 if gene.strand == "+" else gene.end - 500,
        {
            "id": gene.gene_id,
            "strand": gene.strand,
            "row": gene.Index,
            "name": gene.gene_name,
            "type": "gene",
        },
    )
    T[gene.gene_id] = find_tad(gi, tads)
    P[gene.gene_id] = gi
    # add the gene to the GRN for this gene
    G[gene.gene_id].add_node(gi)

# save information in tables, too
P_df = pd.DataFrame(
    {
        "gene_id": [gene_id for gene_id in P.keys()],
        "chr": [prom.chr for prom in P.values()],
        "start": [prom.inf() for prom in P.values()],
        "end": [prom.sup() for prom in P.values()],
    },
)
P_df.to_csv(
    path.join(RESULT_DIR, "promoters.tsv"),
    sep="\t",
    header=True,
    index=False,
)
T_df = pd.DataFrame(
    {
        "gene_id": [gene_id for gene_id in T.keys()],
        "chr": [tad.chr for tad in T.values()],
        "start": [tad.inf() for tad in T.values()],
        "end": [tad.sup() for tad in T.values()],
    },
)
T_df.to_csv(
    path.join(RESULT_DIR, "tads.tsv"),
    sep="\t",
    header=True,
    index=False,
)

logging.info("Adding loops")
L = {}  # loops
for (gene_id, gene) in tqdm(P.items(), total=genes.shape[0]):
    # find loops in this TAD
    tad = T[gene_id]
    tad_loops = loops.loc[
        (
            (loops["chr_x"] == tad.chr)
            & (loops["start_x"] <= tad.sup())
            & (loops["end_x"] >= tad.inf())
        )
        | (
            (loops["chr_y"] == tad.chr)
            & (loops["start_y"] <= tad.sup())
            & (loops["end_y"] >= tad.inf())
        ),
        :,
    ]
    # only include loops where one anchor overlaps the gene gene
    tad_loops = tad_loops.loc[
        # x anchor overlaps the gene (extend to half of loop resolution on either side)
        (
            (tad_loops["start_x"] - LOOP_OFFSET <= gene.sup())
            & (tad_loops["end_x"] + LOOP_OFFSET >= gene.inf())
        )
        # or y anchor overlaps gene
        | (
            (tad_loops["start_y"] - LOOP_OFFSET <= gene.sup())
            & (tad_loops["end_y"] + LOOP_OFFSET >= gene.inf())
        ),
        :,
    ]
    tad_loop_intvls = np.unique(
        [
            Loop(
                l.chr_x,
                l.start_x,
                l.end_x,
                l.chr_y,
                l.start_y,
                l.end_y,
                {
                    "row": l.Index,
                    "id": l.anchor_ID_x,
                },
                {
                    "row": l.Index,
                    "id": l.anchor_ID_y,
                },
                {
                    "row": l.Index,
                    "id": l.loop_ID,
                },
            )
            for l in tad_loops.itertuples()
        ]
    )
    L[gene_id] = tad_loop_intvls.tolist()

# save loops as a table
L_df = pd.concat(
    [
        pd.DataFrame(
            {
                "gene_id": [gene_id for loop in L[gene_id]],
                "chr_x": [loop.left.chr for loop in L[gene_id]],
                "start_x": [loop.left.inf() for loop in L[gene_id]],
                "end_x": [loop.left.sup() for loop in L[gene_id]],
                "anchor_ID_x": [loop.left.data["id"] for loop in L[gene_id]],
                "chr_y": [loop.right.chr for loop in L[gene_id]],
                "start_y": [loop.right.inf() for loop in L[gene_id]],
                "end_y": [loop.right.sup() for loop in L[gene_id]],
                "anchor_ID_y": [loop.right.data["id"] for loop in L[gene_id]],
                "loop_ID": [loop.data["id"] for loop in L[gene_id]],
            }
        )
        for gene_id in L.keys()
    ],
    ignore_index=True,
)
L_df.to_csv(
    path.join(RESULT_DIR, "loops.tsv"),
    sep="\t",
    header=True,
    index=False,
)

logging.info("Linking enhancers to target genes")
E = {}  # Enhancers
for (gene_id, gene) in tqdm(P.items(), total=genes.shape[0]):
    tad = T[gene_id]
    tad_loops = L[gene_id]
    # find all H3K27ac peaks within the TAD
    tad_enhancers = enhancers.loc[
        (enhancers["chr"] == tad.chr)
        & (enhancers["start"] <= tad.sup())
        & (enhancers["end"] >= tad.inf()),
        :,
    ].drop_duplicates()
    # convert to GenomicInterval objects
    tad_enhn_invls = np.unique(
        [
            GenomicInterval(
                e.chr,
                e.start,
                e.end,
                {
                    "row": e.Index,
                    "id": e.Index,
                    "edge": False,
                },
            )
            for e in tad_enhancers.itertuples()
        ]
    )
    # only include enhancers that overlap a loop anchor
    for i, e in enumerate(tad_enhn_invls):
        # add the enhancer to the GRN
        G[gene_id].add_node(e)
        for l in tad_loops:
            # connect the enhancer to the gene if there is a loop connecting them
            # or if they are within the filtering distance from Mustache
            # (the ignored diagonal is 2 x 10 kbp pixels, so 20 kbp)
            if (
                overlapping(e, l.left, LOOP_OFFSET)
                or overlapping(e, l.right, LOOP_OFFSET)
                or overlapping(gene, e, LOOP_DISTANCE_THRESH)
            ):
                # save the edge in the enhancer node itself
                tad_enhn_invls[i].data["edge"] = True
                e = tad_enhn_invls[i]
                G[gene_id].add_edge(gene, e)
    E[gene_id] = tad_enhn_invls.tolist()

E_df = pd.concat(
    [
        pd.DataFrame(
            {
                "gene_id": [gene_id for enh in E[gene_id]],
                "chr": [enh.chr for enh in E[gene_id]],
                "start": [enh.inf() for enh in E[gene_id]],
                "end": [enh.sup() for enh in E[gene_id]],
                "enh_ID": [enh.data["id"] for enh in E[gene_id]],
                "connecting_loop": [enh.data["edge"] for enh in E[gene_id]],
            },
        )
        for gene_id in E.keys()
    ],
    ignore_index=True,
)
# convert to integer types
E_df["start"] = E_df["start"].astype(int)
E_df["end"] = E_df["end"].astype(int)
E_df["enh_ID"] = E_df["enh_ID"].astype(int)
E_df["connecting_loop"] = E_df["connecting_loop"].astype(bool)
E_df.to_csv(
    path.join(RESULT_DIR, "putative-enhancers.tsv"),
    sep="\t",
    header=True,
    index=False,
)

logging.info("Saving graphs as tables")
# convert putative GRNs into a table
# each row is a an edge in the GRN graphs
G_df = genes.merge(right=P_df, on=["gene_id", "chr"], suffixes=["_gene", "_prom"])
G_df = G_df.merge(right=E_df, on=["gene_id", "chr"], suffixes=["", "_enh"])
# only include enhancers that have a loop identified connecting it to the gene
G_df = G_df.loc[G_df["connecting_loop"] == True, :]

# fix column labels and drop unnecessary ones
G_df = G_df.drop(labels=["connecting_loop"], axis=1, inplace=False)
G_df.columns = [
    "chr",
    "start_gene",
    "end_gene",
    "strand",
    "gene_id",
    "gene_name",
    "start_prom",
    "end_prom",
    "start_enh",
    "end_enh",
    "enh_ID",
]

# reorder columns
G_df = G_df.loc[
    :,
    [
        "gene_id",
        "gene_name",
        "chr",
        "start_gene",
        "end_gene",
        "strand",
        "start_prom",
        "end_prom",
        "start_enh",
        "end_enh",
        "enh_ID",
    ],
]

G_df.to_csv(
    path.join(RESULT_DIR, "putative-grns.tsv"),
    sep="\t",
    header=True,
    index=False,
)

logging.info("Calculating GRN Stats")
# for each GRN, count the number gained, shared, and lost loops and enhancers
grn_info = G_df.groupby(["gene_id", "gene_name"]).size().reset_index(name="n_enhancers")

# save count of loops and enhancers in each GRN
grn_info.to_csv(
    path.join(RESULT_DIR, "grn-stats.tsv"),
    sep="\t",
    header=True,
    index=False,
)

# ==============================================================================
# Save data
# ==============================================================================
logging.info("Exporting GRNs")

# save GRNs with pickle
pickle.dump(G, open(path.join(RESULT_DIR, "Graphs", "grns.p"), "wb"))
pickle.dump(T, open(path.join(RESULT_DIR, "Graphs", "tads.p"), "wb"))
pickle.dump(P, open(path.join(RESULT_DIR, "Graphs", "promoters.p"), "wb"))
pickle.dump(E, open(path.join(RESULT_DIR, "Graphs", "enhancers.p"), "wb"))
pickle.dump(L, open(path.join(RESULT_DIR, "Graphs", "loops.p"), "wb"))
