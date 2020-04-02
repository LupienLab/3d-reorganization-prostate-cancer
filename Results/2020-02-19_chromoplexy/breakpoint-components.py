"""
breakpoint-bridging
==========

Take Breakfinder results and link breakpoints together by their positions
"""

import os.path as path
import numpy as np
import pandas as pd
import networkx as nx
from tqdm import tqdm
import pickle
import pprint
from genomic_interval import GenomicInterval, overlapping, find_tad
from copy import deepcopy

# ==============================================================================
# Constants
# ==============================================================================
BREAK_DIR = path.join("..", "2019-07-24_breakfinder", "Breakpoints", "Default")
TAD_DIR = path.join(
    "..", "2020-01-15_TAD-aggregation", "resolved-TADs", "separated-TADs"
)
# distance tolerance for comparisons
TOL = 50000

# window size to check TADs for
W = 3


# ==============================================================================
# Functions
# ==============================================================================
def equivalent_tad(
    bp_i: GenomicInterval,
    bp_j: GenomicInterval,
    tads_i: pd.DataFrame,
    tads_j: pd.DataFrame,
):
    """
    Find if two breakpoints are located in equivalent TADs in their respective samples

    Parameters
    ----------
    bp_i: GenomicInterval
        First breakpoint
    bp_j: GenomicInterval
        Second breakpoint
    tads_i: pd.DataFrame
        TADs from the sample corresponding to `bp_i`
    tads_j : pd.DataFrame
        TADs from the sample corresponding to `bp_j`
    """
    # if the breakpoints are not on the same chromosome, they can't have equivalent TADs
    if bp_i.chr != bp_j.chr:
        return False
    parent_tads_i = find_tad(bp_i, tads_i)
    parent_tads_j = find_tad(bp_j, tads_j)
    # convert to contiguous GenomicInterval to use `overlapping` function
    locus_i = GenomicInterval(
        bp_i.chr, parent_tads_i["start"].min(), parent_tads_i["end"].max()
    )
    locus_j = GenomicInterval(
        bp_j.chr, parent_tads_j["start"].min(), parent_tads_j["end"].max()
    )
    return overlapping(locus_i, locus_j, 0)


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

# load breakpoints and concatenate tables
print("Reading breakpoints")
breakpoints = pd.concat(
    [
        pd.read_csv(
            path.join(BREAK_DIR, s + ".breaks.sorted.manually-resolved.tsv"),
            sep="\t",
            index_col=False,
            names=[
                "chr_x",
                "start_x",
                "end_x",
                "chr_y",
                "start_y",
                "end_y",
                "name",
                "score",
                "strand_x",
                "strand_y",
                "resolution",
                "annotation",
                "notes",
            ],
        )
        for s in SAMPLES
    ],
    keys=SAMPLES,
)

# remove artefacts
breakpoints = breakpoints.loc[breakpoints["annotation"] != "ARTEFACT", :]

# load TADs
print("Reading TADs")
tads = {
    s: pd.read_csv(
        path.join(TAD_DIR, s + ".40000bp.w_" + str(W) + ".domains.bed"),
        sep="\t",
        names=["chr", "start", "end", "lower_persistence", "upper_persistence",],
    )
    for s in SAMPLES
}

# aggregate all TADs together for single search across all samples
all_tads = pd.concat([t for t in tads.values()], keys=SAMPLES)
# move "keys" index to be a column
all_tads = all_tads.reset_index(level=0).rename(columns={"level_0": "Sample"})

# ==============================================================================
# Analysis
# ==============================================================================
# create graph
print("Creating graphs of breakpoints")
G_all = nx.Graph()
G_sample = {s: nx.Graph() for s in SAMPLES}

# add each detected breakpoint as a node to the graph
for s in tqdm(SAMPLES, unit="sample"):
    patient_bps = breakpoints.loc[s]
    for bp in patient_bps.itertuples():
        # create hashable objects to store in each node
        # in this case, and interval
        intvls = [
            GenomicInterval(
                bp.chr_x,
                bp.start_x,
                bp.end_x,
                {"sample": s, "notes": bp.notes, "strand": bp.strand_x, "row": bp.Index},
            ),
            GenomicInterval(
                bp.chr_y,
                bp.start_y,
                bp.end_y,
                {"sample": s, "notes": bp.notes, "strand": bp.strand_y, "row": bp.Index},
            ),
        ]
        # create nodes in the graph
        for i in intvls:
            G_sample[s].add_node(i)
        # link these two nodes since they are linked breakpoints
        G_sample[s].add_edge(intvls[0], intvls[1], annotation=bp.annotation)
print([len(G_sample[s]) for s in SAMPLES])

# breakpoints were detected in pairs, so merge similar ones together
print("Merging redundant breakpoint calls")
pp = pprint.PrettyPrinter(indent=2)
bps_to_merge = []
bp_sets_to_merge = []
for s in tqdm(SAMPLES):
    # redundant breakpoints are merged by calculating the quotient graph, where equivalence is defined by
    # breakpoints that are within 100 kbp of each other and not two ends of a pair (i.e. no edge between them)
    for n in G_sample[s]:
        for m in G_sample[s]:
            if (n == m) or G_sample[s].has_edge(n, m):
                continue
            if overlapping(n, m, TOL // 2):
                n_novel = n not in bps_to_merge
                m_novel = m not in bps_to_merge
                if n_novel and m_novel:
                    bps_to_merge.append(n)
                    bps_to_merge.append(m)
                    bp_sets_to_merge.append(set([n, m]))
                elif m_novel:
                    bps_to_merge.append(m)
                    for stm in bp_sets_to_merge:
                        if n in stm:
                            stm.add(m)
                elif n_novel:
                    bps_to_merge.append(n)
                    for stm in bp_sets_to_merge:
                        if m in stm:
                            stm.add(m)
                #else:
                    #print("Not merging")
                    #continue
pp.pprint([[(i.data["sample"], i.data["row"]) for i in stm] for stm in bp_sets_to_merge])

## manually inspect breakpoints that are in proximity, to see which calls truly need to be merged
# confirmed_bp_sets_to_merge = []
# for stm in bp_sets_to_merge:
#     pp.pprint([(i.data["sample"], i.data["row"], i) for i in stm])
#     print("Merge breakpoints? [y]es to all, [n]o, or [s]ome: ")
#     while True:
#         i = input()
#         if i == "y":
#             confirmed_bp_sets_to_merge.append(stm)
#             break
#         elif i == "n":
#             break
#         else:
#             print("Enter a comma-separated list of breakpoints to merge")
#             stm_list = list(stm)
#             for j, bp in enumerate(stm_list):
#                 print("[" + str(j) + "]" + str(bp))
#             selection = input()
#             selection_ints = [int(j) for j in selection.split(",")]
#             confirmed_bp_sets_to_merge.append(
#                 set([stm_list[j] for j in selection_ints])
#             )
#             break
# pickle.dump(confirmed_bp_sets_to_merge, open("confirmed-breakpoints-to-merge.p", "wb"))

confirmed_bp_sets_to_merge = pickle.load(open("confirmed-breakpoints-to-merge.p", "rb"))


for stm in bp_sets_to_merge:
    G_merged = nx.quotient_graph(
        G=G_sample[s],
        partition=lambda u, v: overlapping(u, v, TOL // 2) and not G_sample[s].has_edge(u, v),
        edge_data=lambda b, c: {
            "annotation": np.unique(
                [G_sample[s].edges[u, v]["annotation"] for u in b for v in c if G_sample[s].has_edge(u, v)]
            )[0] # there should only be one unique SV type (might need to come back to this if it errors)
        }
    )
    # relabel nodes such that the GenomicInterval is the union of all the merged GenomicIntervals
    # each node in G_merged is a frozenset of all the merged nodes
    G_sample[s] = nx.relabel_nodes(
        G_merged,
        {
            nodes: GenomicInterval(
                # they all share the same chromosome, so just get the first value instead of iterating
                np.unique([i.chr for i in nodes])[0],
                # get smallest position
                np.min([i.inf() for i in nodes]),
                # get largest position
                np.max([i.sup() for i in nodes]),
                # choose what information to append to the node as data
                {"sample": s, "notes": [i.data["notes"] for i in nodes], "strand": [i.data["strand"] for i in nodes]}
            ) for nodes in G_merged
        }
    )

print([len(G_sample[s]) for s in SAMPLES])
for s in SAMPLES:
    print(s)
    counts = np.unique([len(cc) for cc in nx.connected_components(G_sample[s])], return_counts=True)
    print(counts[0])
    print(counts[1])

print("Connecting proximal breakpoints")
# create compiled list of all breakpoints for easier parsing
uncoupled_breakpoints = pd.DataFrame(columns=["chr", "start", "end", "mutated_in"])
# connect 2 nodes if their intervals are within 100 kbp of each other
for i, n in tqdm(enumerate(G_all), total=len(G_all)):
    for m in G_all:
        if n == m:
            continue
        # connect these nodes if the identified loci are within 100 kbp of each other
        # the breakpoints are not from the same sample, so this site is recurrent
        # (any nearby points from the same sample have already been merged)
        if overlapping(n, m, TOL // 2):
            G_all.add_edge(n, m, annotation="recurrent")
        # connect these nodes if the identified loci are within the equivalent TADs from their respective samples
        # (I know this isn't the most efficient way to do this, but given the number of breakpoints and samples
        # it's not that much of a concern)
        if equivalent_tad(n, m, tads[n.data["sample"]], tads[m.data["sample"]]):
            G_all.add_edge(n, m, annotation="equivalent-TAD")
    # add index to this node for later use
    n.data["index"] = i
    # save for later (i will be the same as the index in the DataFrame)
    uncoupled_breakpoints = uncoupled_breakpoints.append(
        {"chr": n.chr, "start": n.inf(), "end": n.sup(), "mutated_in": n.data["sample"]},
        ignore_index=True,
        sort=False
    )

# add connectivity of related breakpoints
coupled_tests = pd.DataFrame(columns=["breakpoint_indices", "mutated_in", "n_mut", "n_nonmut"])
for bp in tqdm(G_all, total=len(G_all)):
    # find samples where this, or a nearby, breakpoint occurs
    nbrs = [
        n
        for n, v in G_all[bp].items()
        if v["annotation"] in ["nearby", "recurrent", "equivalent-TAD"]
    ]
    connected_bps = [bp] + nbrs
    indices = ",".join(np.unique([str(n.data["index"]) for n in connected_bps]))
    mut_samples = ",".join(np.unique([n.data["sample"] for n in connected_bps]))
    n_mut = len(np.unique([n.data["sample"] for n in connected_bps]))
    coupled_tests = coupled_tests.append(
        {"breakpoint_indices": indices, "mutated_in": mut_samples, "n_mut": n_mut, "n_nonmut": len(SAMPLES) - n_mut},
        ignore_index=True,
        sort=False
    )
# only keep unique rows
coupled_tests.drop_duplicates(inplace=True, ignore_index=True)

for s in SAMPLES:
    for n in tqdm(G_sample[s]):
        for m in G_sample[s]:
            if n == m:
                continue
            # connect these nodes if the identified loci are within 100 kbp of each other
            if overlapping(n, m, TOL // 2):
                # if the two breakpoints are from the same sample, they're nearby
                G_sample[s].add_edge(n, m, annotation="nearby")
            # connect these nodes if the identified loci at within the same TAD
            if equivalent_tad(n, m, tads[s], tads[s]):
                G_all.add_edge(n, m, annotation="equivalent-TAD")

# ==============================================================================
# Save data
# ==============================================================================
print("Saving graphs")
pickle.dump(G_all, open("breakpoints.all-samples.p", "wb"))
pickle.dump(G_sample, open("breakpoints.per-sample.p", "wb"))
uncoupled_breakpoints.to_csv(
    "sv-breakpoints.tsv",
    sep="\t",
    index_label="breakpoint_index"
)
coupled_tests.to_csv(
    "sv-disruption-tests.tsv",
    sep="\t",
    index_label="test_index"
)

# export to GraphML format
nx.write_graphml(G_all, "Graphs/breakpoints.all-samples.xml")
for s in samples:
    nx.write_graphml(G_sample[s], "Graphs/breakpoints." + s + ".xml")
print("Done")
