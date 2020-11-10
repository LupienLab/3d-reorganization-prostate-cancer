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

# promoter region offset (upstream, downstream)
PROM_OFFSET = (1500, 500)

# get chromosome sizes
hg38 = nc.get_chrominfo("hg38")
CHROM_SIZES = hg38.chrom_lengths

# ==============================================================================
# Functions
# ==============================================================================
def grn_stats(grn_loops) -> Dict[str, int]:
    """
    Count the number of enhancers and loops in each GRN and how they are shared between T2E+/-

    # Parameters
    grn_loops: Dict[str, List[Loop]]
        Loops in the GRN
    grn_enhns: Dict[str, List[GenomicInterval]]
        Enhancers in the GRN
    """
    data = {
        "loops_gained": len(
            [1 for l in grn_loops if l.data["condition"] == "tumour-specific"]
        ),
        "loops_shared": len([1 for l in grn_loops if l.data["condition"] == "shared"]),
        "loops_lost": len(
            [1 for l in grn_loops if l.data["condition"] == "benign-specific"]
        ),
    }
    return data


# ==============================================================================
# Data
# ==============================================================================
logging.info("Loading data")

# load pre-processed promoters
promoters = pd.read_csv("overlapping.promoters.tsv", sep="\t", header=[0],)

# load catalogue of loop calls from all 17 samples
loops = pd.read_csv("overlapping.loops.tsv", sep="\t", header=[0],)

# load sample metadata
metadata = pd.read_csv("config.tsv", sep="\t", header=[0],)
SAMPLES = metadata.loc[metadata.Include == "Yes", "SampleID"].tolist()

# load TAD calls
tads = pd.concat(
    [
        pd.read_csv(
            path.join(
                "..",
                "2020-08-29_TADs-downsampled",
                "Aggregated-TADs",
                "separated-TADs",
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
tads.rename({"level_0": "SampleID"}, axis=1, inplace=True)
tads.drop(labels=["level_1"], axis=1, inplace=True)
# sort by chr, start, end, SampleID
tads.sort_values(
    by=["chr", "start", "end", "SampleID"], inplace=True, ignore_index=True
)


# ==============================================================================
# Analysis
# ==============================================================================
logging.info("Creating GRN")

# create placeholder for the GRN
G = {gid: nx.Graph() for gid in promoters["gene_id"]}

T = {}  # TADs
P = {}  # Promoters
for gene in tqdm(promoters.itertuples(), total=promoters.shape[0]):
    gi = GenomicInterval(gene.chr, gene.start, gene.end)
    T[gene.gene_id] = find_tad(gi, tads)
    if gene.strand == "+":
        prom = GenomicInterval(
            gene.chr,
            max(0, gene.start - PROM_OFFSET[0]),
            min(CHROM_SIZES[gene.chr], gene.start + PROM_OFFSET[1]),
            {
                "row": gene.Index,
                "id": gene.gene_id,
                "name": gene.gene_name,
                "type": "promoter",
            },
        )
    else:
        prom = GenomicInterval(
            gene.chr,
            max(0, gene.end - PROM_OFFSET[1]),
            min(CHROM_SIZES[gene.chr], gene.end + PROM_OFFSET[0]),
            {
                "row": gene.Index,
                "id": gene.gene_id,
                "name": gene.gene_name,
                "type": "promoter",
            },
        )
    P[gene.gene_id] = prom
    # add the promoter to the GRN for this gene
    G[gene.gene_id].add_node(prom)


L = {}  # loops
for (gene_id, prom) in tqdm(P.items(), total=promoters.shape[0]):
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
    # only include loops where one anchor overlaps the gene promoter
    tad_loops = tad_loops.loc[
        # x anchor overlaps the promoter (extend to half of loop resolution on either side)
        (
            (tad_loops["start_x"] - 5000 <= prom.sup())
            & (tad_loops["end_x"] + 5000 >= prom.inf())
        )
        # or y anchor overlaps promoter
        | (
            (tad_loops["start_y"] - 5000 <= prom.sup())
            & (tad_loops["end_y"] + 5000 >= prom.inf())
        ),
        :,
    ]
    tad_loop_intvls = [
        Loop(
            l.chr_x,
            l.start_x,
            l.end_x,
            l.chr_y,
            l.start_y,
            l.end_y,
            {"row": l.Index, "id": l.anchor_ID_x,},
            {"row": l.Index, "id": l.anchor_ID_y,},
            {"row": l.Index, "id": l.loop_ID, "condition": l.loop_type},
        )
        for l in tad_loops.itertuples()
    ]
    L[gene_id] = tad_loop_intvls


logging.info("Calculating GRN Stats")
# for each GRN, count the number gained, shared, and lost loops and enhancers
grn_info = {gene_id: grn_stats(L[gene_id]) for gene_id in G.keys()}
# need to ensure that the gene_ids are always listed in the same order
grn_cre_stats = pd.DataFrame(
    {
        "gene_id": [gene_id for gene_id in promoters["gene_id"]],
        "loops_gained": [
            grn_info[gene_id]["loops_gained"] for gene_id in promoters["gene_id"]
        ],
        "loops_shared": [
            grn_info[gene_id]["loops_shared"] for gene_id in promoters["gene_id"]
        ],
        "loops_lost": [
            grn_info[gene_id]["loops_lost"] for gene_id in promoters["gene_id"]
        ],
    }
)


# ==============================================================================
# Save data
# ==============================================================================
logging.info("Exporting GRN")

# save GRNs with pickle
pickle.dump(G, open(path.join("Graphs", "grns.p"), "wb"))
pickle.dump(T, open(path.join("Graphs", "tads.p"), "wb"))
pickle.dump(P, open(path.join("Graphs", "promoters.p"), "wb"))
pickle.dump(L, open(path.join("Graphs", "loops.p"), "wb"))

# save count of loops and enhancers in each GRN
grn_cre_stats.to_csv("Graphs/grn-stats.tsv", sep="\t", index=False)
