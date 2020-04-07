"""
breakpoint-bridging
==========

Take Breakfinder results and link breakpoints together by their positions
"""

from typing import List, Set, Tuple
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
GRAPH_DIR = "Graphs"
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


def merge_nodes(G: nx.Graph, nodes: List[GenomicInterval]) -> nx.Graph:
    # relabel nodes such that the GenomicInterval is the union of all the merged GenomicIntervals
    merged_node = GenomicInterval(
        chrom=nodes[0].chr,
        # get smallest position
        start=np.min([i.inf() for i in nodes]),
        # get largest position
        end=np.max([i.sup() for i in nodes]),
        # choose what information to append to the node as data
        data={
            "sample": nodes[0].data["sample"],
            "notes": [i.data["notes"] for i in nodes],
            "strand": [i.data["strand"] for i in nodes]
        }
    )
    G.add_node(merged_node)
    # iterate over all edges, and copy them to merged_node the edge connects to one of the nodes
    edges_to_add = {}
    for n1, n2, data in G.edges(data=True):
        n1_being_merged = n1 in nodes
        n2_being_merged = n2 in nodes
        if n1_being_merged and n2_being_merged:
            # don't create self loops, since this is the same breakpoint
            continue
        # can't modify G while iterating over it, so the old edges and nodes have to be
        # modified afterwards
        elif n1 in nodes:
            edges_to_add[n2] = data
        elif n2 in nodes:
            edges_to_add[n1] = data
    
    # add edges
    for node, data in edges_to_add.items():
        G.add_edge(merged_node, node, **data)
    # remove old nodes
    for n in nodes:
        # removing a node also removes all edges connected to that node
        G.remove_node(n)
    return G


def compatible_sv_pair(a: str, b: str) -> bool:
    """
    Determine whether a pair of structural variant types are compatible
    """
    always_compatible = ["BND", "UNKNOWN"]
    if a in always_compatible:
        return True
    elif b in always_compatible:
        return True
    else:
        return a == b


def compatible_sv_types(type_a: Set[str], type_b: Set[str]) -> Tuple[bool, Set[str]]:
    """
    List all compatible structural variant types between two sets
    """
    # find whether all pairs of SV types are compatible, return if they are not
    for a in type_a:
        for b in type_b:
            if not compatible_sv_pair(a, b):
                return (False, set([]))
    # if no incompatible pair has been found, return True
    return (True, type_a.union(type_b))
    

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
# 1. Create graphs
# --------------------------------------
print("Creating graphs of breakpoints")
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
        for intvl in intvls:
            G_sample[s].add_node(intvl)
        # link these two nodes since they are linked breakpoints
        G_sample[s].add_edge(intvls[0], intvls[1], annotation=bp.annotation)

# 2. Merge redundant breakpoints
# --------------------------------------
# hic_breakfinder identifies breakpoint pairs, so the same breakpoint can be called multiple times
# if it interacts with multiple other breakpoints
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

#pp.pprint([[(i.data["sample"], i.data["row"]) for i in stm] for stm in bp_sets_to_merge])

# 2. a) manually inspect breakpoints that are in proximity, to see which calls truly need to be merged
# --------------------------------------

# this commented out section is run once, manually, to inspect each set of nearby breakpoints
# this list is then saved to "confirmed-breakpoints-to-merge.p"

# confirmed_bp_sets_to_merge = []
# for stm in bp_sets_to_merge:
#     pp.pprint([(i.data["sample"], i.data["row"], i) for i in stm])
#     print("Merge breakpoints? [y]es to all, [n]o, or [s]ome: ")
#     while True:
#         user_input = input()
#         if user_input == "y":
#             confirmed_bp_sets_to_merge.append(stm)
#             break
#         elif user_input == "n":
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
# pickle.dump(
#     [G_sample, confirmed_bp_sets_to_merge],
#     open(path.join(GRAPH_DIR, "breakpoints.per-sample.with-multiplicity.p"), "wb")
# )

# after the manual confirmation of these breakpoints, the redundant breakpoints are loaded
# note that these objects have to be saved and loaded together, because of how networkx stores nodes in its graphs
G_sample_copy, confirmed_bp_sets_to_merge_copy = pickle.load(open(
    path.join(GRAPH_DIR, "breakpoints.per-sample.with-multiplicity.p"),
    "rb"
))
G_sample = G_sample_copy
confirmed_bp_sets_to_merge = confirmed_bp_sets_to_merge_copy

# 2. b) contract the redundant breakpoints in the graphs
# --------------------------------------
for stm in confirmed_bp_sets_to_merge:
    nodes = list(stm)
    s = nodes[0].data["sample"]
    G_sample[s] = merge_nodes(G_sample[s], nodes)

# add a unique ID to the breakpoint in the dataset for future reference
bp_counter = 0
for s in SAMPLES:
    for n in G_sample[s]:
        n.data["breakpoint_ID"] = bp_counter
        bp_counter += 1

# add an ID for the component to which each breakpoint belongs (each connected component is a different SV event)
for s in SAMPLES:
    for component_counter, cc in enumerate(nx.connected_components(G_sample[s])):
        for n in cc:
            n.data["component_ID"] = component_counter

# 2. c) save merged breakpoints in various formats
# --------------------------------------
# save merged breakpoints in pickle file
pickle.dump(G_sample, open(path.join(GRAPH_DIR, "breakpoints.per-sample.merged-breakpoints.p"), "wb"))

# create compiled list of all breakpoints for easier parsing and saving in a table
uncoupled_breakpoints = pd.DataFrame(columns=["chr", "start", "end", "SampleID", "breakpoint_ID", "component_ID"])
for s in SAMPLES:
    for n in G_sample[s]:
        uncoupled_breakpoints = uncoupled_breakpoints.append(
            {
                "chr": n.chr,
                "start": n.inf(),
                "end": n.sup(),
                "SampleID": n.data["sample"],
                "breakpoint_ID": n.data["breakpoint_ID"],
                "component_ID": n.data["component_ID"],
            },
            ignore_index=True,
            sort=False
        )

# save in a table
uncoupled_breakpoints.to_csv(
    path.join(GRAPH_DIR, "sv-breakpoints.tsv"),
    sep="\t",
    index=False
)

# produce global list of merged breakpoints
paired_breakpoints = pd.DataFrame(
    columns=[
        "chr_x", "start_x", "end_x",
        "chr_y", "start_y", "end_y",
        "breakpoint_ID_x", "breakpoint_ID_y",
        "component_ID_x", "component_ID_y",
        "SampleID", "sv_type"
    ]
)
for s in SAMPLES:
    for n1, n2, data in G_sample[s].edges(data=True):
        # sort the nodes prior to listing, so the lesser of the two breakpoints is always in the first column
        nodes = sorted([n1, n2])
        paired_breakpoints = paired_breakpoints.append(
            {
                "chr_x": nodes[0].chr,
                "start_x": nodes[0].inf(),
                "end_x": nodes[0].sup(),
                "chr_y": nodes[1].chr,
                "start_y": nodes[1].inf(),
                "end_y": nodes[1].sup(),
                "breakpoint_ID_x": nodes[0].data["breakpoint_ID"],
                "breakpoint_ID_y": nodes[1].data["breakpoint_ID"],
                "component_ID_x": nodes[0].data["component_ID"],
                "component_ID_y": nodes[1].data["component_ID"],
                "SampleID": s,
                "sv_type": data["annotation"]
            },
            ignore_index = True,
            sort=False
        )

# save as table
paired_breakpoints.to_csv(
    path.join(GRAPH_DIR, "sv-breakpoints.paired.tsv"),
    sep="\t",
    index=False
)

# 3. Connecting points in graph containing all samples based on their TADs and locations
# --------------------------------------
print("Connecting proximal breakpoints")
# create graph containing all samples
G_all = nx.Graph()
for s in SAMPLES:
    for n1, n2, data in G_sample[s].edges(data=True):
        G_all.add_edge(n1, n2, **data)

# determine which breakpoints need to be grouped together for hypothesis testing
test_ID_status = {}
test_IDs_to_merge = []
test_IDs_to_never_merge = []
test_ID_SV_types = {}
for i, n in tqdm(enumerate(G_all), total=len(G_all)):
    # find all SV types that this breakpoint is involved in
    n_sv_types = set([e["annotation"] for e in G_all[n].values()])
    # determine whether this breakpoint has already been assigned to a testing group
    n_previously_assigned_test = "test_ID" in n.data
    if not n_previously_assigned_test:
        n.data["test_ID"] = i
        test_ID_status[i] = True
        test_ID_SV_types[i] = set([e["annotation"] for e in G_all[n].values()])
    for j, m in enumerate(G_all):
        # iteration is always in the same order, skip any node that has already been assessed with this one
        if j <= i:
            continue
        # connect these nodes if the identified loci are within the equivalent TADs from their respective samples
        # (I know this isn't the most efficient way to do this, but given the number of breakpoints and samples
        # it's not that much of a concern)
        if overlapping(m, n, TOL // 2) or equivalent_tad(n, m, tads[n.data["sample"]], tads[m.data["sample"]]):
            # find all SV types that this breakpoint is involved in
            m_sv_types = set([e["annotation"] for e in G_all[m].values()])
            # determine whether this breakpoint has already been assigned to a testing group
            m_previously_assigned_test = "test_ID" in m.data
            # determine how to group these breakpoints for hypothesis testing
            is_compatible, compatible_types = compatible_sv_types(n_sv_types, m_sv_types)
            if is_compatible:
                # if there is a compatible SV types between them
                # (e.g. they both have DEL in the same TAD)
                # group these breakpoints together in the same testing group and test them
                if m_previously_assigned_test:
                    # if both nodes have already been assigned test groups, be sure to merge them
                    if n.data["test_ID"] != m.data["test_ID"]:
                        test_IDs_to_merge.append([n.data["test_ID"], m.data["test_ID"]])
                else:
                    # if m has not been assigned a test group yet, assign it to n's test group
                    m.data["test_ID"] = n.data["test_ID"]
                    if n.data["test_ID"] not in test_ID_status:
                        test_ID_status[n.data["test_ID"]] = True
                if n.data["test_ID"] in test_ID_SV_types:
                    test_ID_SV_types[n.data["test_ID"]] = test_ID_SV_types[n.data["test_ID"]].union(compatible_types)
                else:
                    test_ID_SV_types[n.data["test_ID"]] = compatible_types
            elif n.data["sample"] != m.data["sample"]:
                # if different events from different samples
                # group separately and test separately
                if not m_previously_assigned_test:
                    m.data["test_ID"] = j
                    test_ID_status[j] = True
                    test_ID_SV_types[j] = set([e["annotation"] for e in G_all[m].values()])
                test_IDs_to_never_merge.append([n.data["test_ID"], m.data["test_ID"]])
            else:
                # group these breakpoints together in the same testing group but don't test them
                # can't disambiguate the effect of one event vs the other at this region
                if m_previously_assigned_test:
                    if n.data["test_ID"] != m.data["test_ID"]:
                        test_IDs_to_merge.append([n.data["test_ID"], m.data["test_ID"]])
                else:
                    m.data["test_ID"] = n.data["test_ID"]
                # this sets this entire testing group to False, if any of its members are ambiguous
                test_ID_status[n.data["test_ID"]] = False

# merge testing groups that were initially given separate IDs but were later discovered to be part of the same group
# if the pair of test groups are never to be merged (i.e. are ambiguous), be sure to separate the samples for background tests
ttm_graph = nx.Graph()
ttnm_graph = nx.Graph()
ttm_new_IDs = {}
for ttm in test_IDs_to_merge:
    ttm_graph.add_edge(ttm[0], ttm[1])
for ttnm in test_IDs_to_never_merge:
    ttnm_graph.add_edge(ttnm[0], ttnm[1])

# reduce each test_ID to the smallest in that component
for cc in nx.connected_components(ttm_graph):
    min_test_ID = min(n for n in cc)
    for n in cc:
        ttm_new_IDs[n] = min_test_ID
        # take the union of SV types across these similar test groups
        if n in test_ID_SV_types:
            test_ID_SV_types[min_test_ID] = test_ID_SV_types[min_test_ID].union(test_ID_SV_types[n])

# apply this reduction to the nodes that should not be merged (i.e. distinct SV events proximal to each other)
ttnm_graph = nx.relabel_nodes(ttnm_graph, ttm_new_IDs)
# change testing status for any test group that cannot be merged with itself
for e in ttnm_graph.edges():
    if e[0] == e[1]:
        # this may already be the case, but this ensures it
        test_ID_status[e[0]] = False

# finalize unique tests
reduced_test_IDs = set([])
for n in G_all:
    if n.data["test_ID"] in ttm_new_IDs:
        # if one group is ambiguous and not to be tested, then the aggregated group should not be either
        test_ID_status[ttm_new_IDs[n.data["test_ID"]]] = test_ID_status[ttm_new_IDs[n.data["test_ID"]]] and test_ID_status[n.data["test_ID"]]
        # re-label the test_IDs for any test in test_IDs_to_merge
        n.data["test_ID"] = ttm_new_IDs[n.data["test_ID"]]
    reduced_test_IDs.add(n.data["test_ID"])

# remove unnecessary test IDs
test_ID_status = {k:v for k, v in test_ID_status.items() if k in reduced_test_IDs}

# save graph to pickle object
pickle.dump(G_all, open(path.join(GRAPH_DIR, "breakpoints.all-samples.p"), "wb"))

# create table of test IDs and metadata
coupled_tests = pd.DataFrame(
    columns=[
        "test_ID",
        "breakpoint_IDs",
        "mut_samples",
        "nonmut_samples",
        "n_mut",
        "n_nonmut",
        "SV_type",
        "testing",
    ]
)
for t_id in test_ID_status:
    # get IDs of all breakpoints in this test group
    b_ids = set([n.data["breakpoint_ID"] for n in G_all if n.data["test_ID"] == t_id])
    # get all the samples involved in this test group
    mut_samples = set([n.data["sample"] for n in G_all if n.data["test_ID"] == t_id])
    # get all the non-mutated samples that will be used as the null distribution
    nonmut_samples = set(SAMPLES) - mut_samples
    # be sure to remove any samples in the tests that are incompatible with this one
    for e in ttnm_graph.edges():
        incompatible_samples = set([])
        if t_id == e[0]:
            incompatible_samples = set([n.data["sample"] for n in G_all if n.data["test_ID"] == e[1]])
        elif t_id == e[1]:
            incompatible_samples = set([n.data["sample"] for n in G_all if n.data["test_ID"] == e[0]])
        nonmut_samples = nonmut_samples - incompatible_samples
    coupled_tests = coupled_tests.append(
        {
            "test_ID": t_id,
            "breakpoint_IDs": ",".join([str(i) for i in b_ids]),
            "mut_samples": ",".join(sorted(mut_samples)),
            "nonmut_samples": ",".join(sorted(nonmut_samples)),
            "n_mut": len(mut_samples),
            "n_nonmut": len(nonmut_samples),
            "SV_type": ",".join(sorted(test_ID_SV_types[t_id])),
            "testing": test_ID_status[t_id]
        },
        ignore_index=True,
        sort=False
    )

# save coupled tests
coupled_tests.to_csv(
    path.join(GRAPH_DIR, "sv-disruption-tests.tsv"),
    sep="\t",
    index=False
)

# ==============================================================================
# Save data
# ==============================================================================
print("Exporting graphs to GraphML")
# export to GraphML format
nx.write_graphml(G_all, "Graphs/breakpoints.all-samples.xml")
for s in SAMPLES:
    nx.write_graphml(G_sample[s], "Graphs/breakpoints." + s + ".xml")
print("Done")
