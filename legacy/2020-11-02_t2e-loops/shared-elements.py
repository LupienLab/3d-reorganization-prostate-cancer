# ==============================================================================
# Meta
# ==============================================================================
# shared-elements
# --------------------------------------
# Description: Test if differentially expressed genes from the same SV are more likely to share enhancers, or not
# Author: James Hawley

import os.path as path
import numpy as np
import pandas as pd
import pickle
import networkx as nx
from tqdm import tqdm

from genomic_interval import GenomicInterval, overlapping

# ==============================================================================
# Functions
# ==============================================================================
def count_shared_enhancers(gene_id_1, gene_id_2):
    # get the nodes in the GRNs that are the genes themselves, not the enhancers
    gene_1_body = [el for el in G[gene_id_1].nodes() if el.data["type"] == "gene"][0]
    gene_2_body = [el for el in G[gene_id_2].nodes() if el.data["type"] == "gene"][0]
    # get the enhancers for the two genes that connect directly to the gene
    g1_enhns = G[gene_id_1][gene_1_body]
    g2_enhns = G[gene_id_2][gene_2_body]
    # count how many are shared
    # it is sufficient to count the intersection of the enhancer IDs for the two genes
    g1_enhns_ids = set([e.data["id"] for e in g1_enhns])
    g2_enhns_ids = set([e.data["id"] for e in g2_enhns])
    return len(g1_enhns_ids.intersection(g2_enhns_ids))


# ==============================================================================
# Data
# ==============================================================================
sv_exprs = pd.read_csv(
    "../2020-06-18_sv-disruption-expression/summary-sv-disruption.tsv", sep="\t"
)
G = pickle.load(open("Graphs/grns.p", "rb"))
T = pickle.load(open("Graphs/tads.p", "rb"))

gene_ids = list(G.keys())
gene_names = {event.target_id: event.gene_name for event in sv_exprs.itertuples()}
event_IDs = np.unique(sv_exprs["event_ID"]).tolist()

# ==============================================================================
# Analysis
# ==============================================================================
# 1. Get genes that share a TAD within the same SV
# ------------------------------------------------
linked_sv_genes = {eid: nx.Graph() for eid in event_IDs}

for eid in tqdm(event_IDs):
    # get all genes in that event
    event_genes = sv_exprs.loc[sv_exprs["event_ID"] == eid, ["target_id", "qval"]]
    event_genes["dge"] = event_genes["qval"] < 0.05
    gene_ids = event_genes["target_id"].tolist()
    dge = {r.target_id: r.dge for r in event_genes.itertuples()}
    for i, gi in enumerate(gene_ids):
        if gi not in T.keys():
            continue
        tad_i = T[gi]
        linked_sv_genes[eid].add_node(gi, dge=dge[gi])
        for gj in gene_ids[:i]:
            if gj not in T.keys():
                continue
            tad_j = T[gj]
            if overlapping(tad_i, tad_j):
                linked_sv_genes[eid].add_node(gj)
                linked_sv_genes[eid].add_edge(gi, gj)


# 2. Calculate the number of shared enhancers between any two genes
# -----------------------------------------------------------------
for eid in tqdm(event_IDs):
    for gene1, gene2 in linked_sv_genes[eid].edges():
        shared_enhns = count_shared_enhancers(gene1, gene2)
        # add this number to the edge in the graph
        linked_sv_genes[eid][gene1][gene2]["shared_enhancers"] = shared_enhns


# 3. Calculate the mean number of shared enhancers between two genes that are both differentially expressed or not
# ----------------------------------------------------------------------------------------------------------------
obs_shared_enhancers = pd.DataFrame(
    columns=[
        "event_ID",
        "gene_id_1",
        "gene_name_1",
        "gene_id_2",
        "gene_name_2",
        "both_DGE",
        "N_shared_enhancers",
    ]
)

for eid in tqdm(event_IDs):
    dge_status = nx.get_node_attributes(linked_sv_genes[eid], "dge")
    for gene1, gene2, edata in linked_sv_genes[eid].edges(data=True):
        obs_shared_enhancers = obs_shared_enhancers.append(
            {
                "event_ID": eid,
                "gene_id_1": gene1,
                "gene_name_1": gene_names[gene1],
                "gene_id_2": gene2,
                "gene_name_2": gene_names[gene2],
                "both_DGE": (dge_status[gene1] and dge_status[gene2]),
                "N_shared_enhancers": edata["shared_enhancers"],
            },
            ignore_index=True,
        )


# 4. Calculate median of shared enhancers between the gene pairs
# --------------------------------------------------------------
shared_enhn_count = obs_shared_enhancers["N_shared_enhancers"].tolist()
n_both_dge_comps = obs_shared_enhancers.loc[
    obs_shared_enhancers["both_DGE"] == True, :
].shape[0]

obs_log2fc = np.log2(
    np.mean(
        obs_shared_enhancers.loc[
            obs_shared_enhancers["both_DGE"] == True, "N_shared_enhancers"
        ]
    )
    / np.mean(
        obs_shared_enhancers.loc[
            obs_shared_enhancers["both_DGE"] == False, "N_shared_enhancers"
        ]
    )
)

n_perms = 10000
perms = pd.DataFrame(
    {"Mean_Both_DGE_Enhns": [0] * n_perms, "Mean_Some_NonDGE_Enhns": [0] * n_perms,}
)
# perform permutation
for i in tqdm(range(n_perms)):
    np.random.shuffle(shared_enhn_count)
    perms.iloc[i]["Mean_Both_DGE_Enhns"] = np.median(
        shared_enhn_count[0:n_both_dge_comps]
    )
    perms.iloc[i]["Mean_Some_NonDGE_Enhns"] = np.median(
        shared_enhn_count[n_both_dge_comps:]
    )

perms["log2FC"] = np.log2(
    perms["Mean_Both_DGE_Enhns"] / perms["Mean_Some_NonDGE_Enhns"]
)

pval = perms.loc[perms["log2FC"] > obs_log2fc].shape[0] / n_perms
print(pval)

# ==============================================================================
# Save data
# ==============================================================================
perms.to_csv(
    "shared-enhancers.permutations.tsv", sep="\t", index=True, index_label="Iteration"
)

obs_shared_enhancers.to_csv("shared-enhancers.tsv", sep="\t", index=False)
