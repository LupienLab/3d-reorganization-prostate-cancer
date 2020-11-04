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
def sat_grn(grn: nx.Graph) -> bool:
    """
    Determine if a GRN for a specific gene is satisfactory for the analysis involving genes, enhancers, and expression

    Parameters
    ==========
    grn: Gene regulatory network for a single gene
    """
    loop_conditions = set(
        [loop_data["condition"] for _, _, loop_data in grn.edges(data=True)]
    )
    enhn_conditions = set(
        [e.data["condition"] for e in grn if e.data["type"] == "enhancer"]
    )
    # if all the loops and enhancers are shared, don't test it
    all_shared_loops = loop_conditions == set(["shared"])
    all_shared_enhns = enhn_conditions == set(["shared"])
    if all_shared_loops and all_shared_enhns:
        return False
    elif all_shared_loops:
        if ("T2E-specific" in enhn_conditions) and (
            "nonT2E-specific" in enhn_conditions
        ):
            return False
    elif all_shared_enhns:
        if ("T2E-specific" in loop_conditions) and (
            "nonT2E-specific" in loop_conditions
        ):
            return False
        return False
    # if there are multiple loops with different changes of direction (i.e. one TE2+ and one T2E- loop), skip
    if ("T2E-specific" in loop_conditions) and ("nonT2E-specific" in loop_conditions):
        return False
    # if the direction of change in the loop and the enhancer are opposite, don't test
    for cre1, cre2, loop_data in grn.edges(data=True):
        if cre1.data["type"] == "enhancer":
            enhn_cdn = cre1.data["condition"]
        else:
            enhn_cdn = cre2.data["condition"]
        if set([loop_data["condition"], enhn_cdn]) == set(
            ["T2E-specific", "nonT2E-specific"]
        ):
            return False
    return True


# ==============================================================================
# Data
# ==============================================================================
logging.info("Loading data")

# load pre-processed promoters
promoters = pd.read_csv("overlapping.promoters.tsv", sep="\t", header=[0],)

# load catalogue of H3K27ac peaks and differential test results
enhancers = pd.read_csv("overlapping.enhancers.tsv", sep="\t", header=[0],)

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
for (gene_id, prom), (_, tad) in tqdm(
    zip(P.items(), T.items()), total=promoters.shape[0]
):
    # find loops in this TAD
    tad_loops = loops.loc[
        (loops["chr_x"] == tad.chr)
        & (loops["chr_y"] == tad.chr)
        & (loops["start_x"] <= tad.sup())
        & (loops["start_y"] <= tad.sup())
        & (loops["end_x"] <= tad.inf())
        & (loops["end_y"] <= tad.inf()),
        :,
    ]
    # only include loops where one anchor overlaps the gene promoter
    tad_loops = tad_loops.loc[
        # x anchor overlaps the promoter
        ((tad_loops["start_x"] <= prom.sup()) & (tad_loops["end_x"] <= prom.inf()))
        # or y anchor overlaps promoter
        | ((tad_loops["start_y"] <= prom.sup()) & (tad_loops["end_y"] <= prom.inf())),
        :,
    ]
    # no loops connecting to gene, no GRNs to be made, so skip the rest
    if len(tad_loops) == 0:
        continue
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
            {"row": l.Index, "id": l.loopID, "condition": l.loop_type},
        )
        for l in tad_loops.itertuples()
    ]
    L[gene_id] = tad_loop_intvls


E = {}  # Enhancers
for (gene_id, prom), (_, tad), (_, tad_loops) in tqdm(
    zip(P.items(), T.items(), L.items()), total=promoters.shape[0]
):
    # find all H3K27ac peaks within the TAD, exlcuding promoter regions
    tad_enhancers = enhancers.loc[
        (enhancers.index != prom.data["row"])
        & (enhancers["chr"] == tad.chr)
        & (enhancers["start"] <= tad.sup())
        & (enhancers["end"] <= tad.inf()),
        :,
    ]
    # if no enhancers, no GRN to be made, so skip the rest
    if len(tad_enhancers) == 0:
        continue
    # calculate whether the enhancer is condition-specific or shared
    enhn_conditions = [""] * tad_enhancers.shape[0]
    for i, e in enumerate(tad_enhancers.itertuples()):
        if e.FDR < 0.05:
            if e.Fold > 0:
                enhn_conditions[i] = "T2E-specific"
            else:
                enhn_conditions[i] = "nonT2E-specific"
        else:
            enhn_conditions[i] = "shared"
    # convert to GenomicInterval objects
    tad_peak_intvls = [
        GenomicInterval(
            e.chr,
            e.start,
            e.end,
            {
                "row": e.Index,
                "id": e.Index,
                "name": "",
                "type": "enhancer",
                "sig": e.FDR < 0.05,
                "fc": e.Fold,
                "condition": cdn,
            },
        )
        for p, cdn in zip(tad_enhancers.itertuples(), enhn_conditions)
    ]
    # only include enhancers that overlap a loop anchor
    for e in tad_peak_intvls:
        # this may produce multiple edges between the enhancer and the promoters
        # (although biologically, this shouldn't happen at the enhancers, only the promoters)
        # I will need to check for this in future code
        for l in tad_loops:
            if overlapping(e, l.left) or overlapping(e, l.right):
                G[gene_id].add_node(e)
                # connect the enhancer to the gene promoter
                G[gene_id].add_edge(prom, e, condition=l.data["condition"])
    E[gene_id] = tad_peak_intvls


logging.info("Checking GRN satisfiability")
# only consider genes that have at least 1 loop and at least 1 enhancer
genes_with_loops = set(
    [gene_id for gene_id, grn_loops in L.items() if len(grn_loops) > 0]
)
genes_with_enhancers = set(
    [gene_id for gene_id, grn_enhns in E.items() if len(grn_enhns) > 0]
)
genes_to_consider = genes_with_loops.intersection(genes_with_enhancers)
# remove other genes from the dict of GRNs
G = {gene_id: grn for gene_id, grn in G.items() if gene_id in genes_to_consider}
T = {gene_id: tad for gene_id, tad in T.items() if gene_id in genes_to_consider}
P = {gene_id: prom for gene_id, prom in P.items() if gene_id in genes_to_consider}
E = {gene_id: enhn for gene_id, enhn in E.items() if gene_id in genes_to_consider}
L = {gene_id: loop for gene_id, loop in L.items() if gene_id in genes_to_consider}

grn_sat_table = pd.DataFrame(
    {
        "gene_id": [gene_id for gene_id in G.keys()],
        "Satisfies_Testing_Requirements": [False] * len(G),
    }
)
# check if each gene has a satisfactory GRN
for gene in tqdm(grn_sat_table.itertuples(), total=grn_sat_table.shape[0]):
    # if the GRN is good for unambiguous relationships between genes, loops, and enhancers
    if sat_grn(G[gene.gene_id]):
        grn_sat_table.loc[gene.Index, "Satisfies_Testing_Requirements"] = True

# ==============================================================================
# Save data
# ==============================================================================
logging.info("Exporting GRN")

# save GRNs with pickle
pickle.dump(G, open(path.join("Graphs", "grns.p"), "wb"))
pickle.dump(T, open(path.join("Graphs", "tads.p"), "wb"))
pickle.dump(P, open(path.join("Graphs", "promoters.p"), "wb"))
pickle.dump(E, open(path.join("Graphs", "enhancers.p"), "wb"))
pickle.dump(L, open(path.join("Graphs", "loops.p"), "wb"))

# save the satisfiability table for future reference
grn_sat_table.to_csv(
    "Graphs/grn-satisfiability.tsv", sep="\t",
)

