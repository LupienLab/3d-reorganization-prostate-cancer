"""
compare-sv
==========

Compare SVs identified with WGS + Delly to the SVs identified by Hi-C + hic_breakfinder
"""

import os.path as path
import numpy as np
import pandas as pd
from tqdm import tqdm
import itertools

RES_DIR = path.join("..", "..", "results", "2020-02-20_wgs-hic-sv-comparison")


# ==============================================================================
# Constants
# ==============================================================================
# tolerance for interval overlap
TOL = 50000

# ==============================================================================
# Data
# ==============================================================================
# load metadata
metadata = pd.read_csv(
    path.join("..", "..", "data", "External", "LowC_Samples_Data_Available.tsv"),
    sep="\t",
    header=0,
    index_col=False,
)
metadata = metadata.loc[metadata.Include == "Yes",:]
metadata["SampleID"] = ["PCa" + str(i) for i in metadata["Sample ID"]]
SAMPLES = metadata["SampleID"].tolist()

# load WGS breakpoints
wgs = pd.concat(
    [
        pd.read_csv(
            path.join(
                "..", "..", "data", "External","CPC-GENE",
                "structural-variant-vcfs", s + ".hg38.sorted.bed"
            ),
            sep="\t",
            header=None,
            names=["chr", "start", "end", "SV_ID"],
        ) for s in SAMPLES
    ],
    keys=SAMPLES,
    names=["SampleID", "i"],
)
wgs.reset_index(inplace=True)
wgs.drop(labels="i", axis=1, inplace=True)
wgs["breakpoint_ID"] = (
    wgs["SampleID"].astype(str)
    + "_"
    + wgs["SV_ID"].astype(str)
    + "_"
    + wgs.index.astype(str)
)

# load Hi-C breakpoints
hic = pd.read_csv(
    path.join(
        "..", "..", "results",
        "2020-02-19_chromoplexy", "Graphs",
        "sv-breakpoints.tsv"
    ),
    sep="\t",
    header=0,
    index_col=False,
)

# ==============================================================================
# Analysis
# ==============================================================================
wgs["Mutual"] = False
wgs["Mutual_Breakpoint_IDs"] = ""
hic["Mutual"] = False
hic["Mutual_Breakpoint_IDs"] = ""
for s in SAMPLES:
    # reduce tables to the specific sample for faster queries
    sample_wgs = wgs.loc[wgs.SampleID == s, :]
    sample_hic = hic.loc[hic.SampleID == s,:]
    # for each breakpoint called from the WGS data, try to identify a similar
    # breakpoint(s) from the Hi-C data
    for r in sample_wgs.itertuples():
        similar_hic_breaks = sample_hic.loc[
            (
                (sample_hic.chr == r.chr)
                & (r.start - TOL <= sample_hic.end)
                & (r.end + TOL >= sample_hic.start)
            ),
            :
        ]
        if similar_hic_breaks.shape[0] > 0:
            wgs.loc[r.Index, "Mutual"] = True
            wgs.loc[r.Index, "Mutual_Breakpoint_IDs"] = ",".join([
                str(i) for i in similar_hic_breaks.breakpoint_ID
            ])
    # for each breakpoint called from the Hi-C data, try to identify a similar
    # breakpoint(s) from the WGS data
    for r in sample_hic.itertuples():
        similar_wgs_breaks = sample_wgs.loc[
            (
                (sample_wgs.chr == r.chr)
                & (r.start <= sample_wgs.end + TOL)
                & (r.end >= sample_wgs.start - TOL)
            ),
            :
        ]
        if similar_wgs_breaks.shape[0] > 0:
            hic.loc[r.Index, "Mutual"] = True
            hic.loc[r.Index, "Mutual_Breakpoint_IDs"] = ",".join([
                str(i) for i in similar_wgs_breaks.breakpoint_ID
            ])

# summarize mutual detections across entire cohort
detections = pd.DataFrame({
    "Source": ["WGS", "Hi-C"],
    "Detected_In_Source": [wgs.shape[0], hic.shape[0]],
    "Mutually_Detected": [
        wgs.loc[wgs.Mutual == True].shape[0],
        hic.loc[hic.Mutual == True].shape[0],
    ],
})

# summarize mutual detections per sample
detections_by_sample = pd.DataFrame({
    "SampleID": list(itertools.chain.from_iterable(
        itertools.repeat(s, 2) for s in SAMPLES
    )),
    "Source": ["WGS", "Hi-C"] * len(SAMPLES),
    "Detected_In_Source": list(itertools.chain.from_iterable(
        [[
            wgs.loc[wgs.SampleID == s].shape[0],
            hic.loc[hic.SampleID == s].shape[0]
        ] for s in SAMPLES]
    )),
    "Mutually_Detected": list(itertools.chain.from_iterable(
        [
            [
                wgs.loc[(wgs.SampleID == s) & (wgs.Mutual == True)].shape[0],
                hic.loc[(hic.SampleID == s) & (hic.Mutual == True)].shape[0],
            ] for s in SAMPLES
        ]
    )),
})


# ==============================================================================
# Save data
# ==============================================================================
detections.to_csv(
    path.join(RES_DIR, "detections.all.tsv"),
    sep="\t",
    index=False
)
detections_by_sample.to_csv(
    path.join(RES_DIR, "detections.per-sample.tsv"),
    sep="\t",
    index=False
)
wgs.to_csv(
    path.join(RES_DIR, "detections.wgs.tsv"),
    sep="\t",
    index=False
)
hic.to_csv(
    path.join(RES_DIR, "detections.hic.tsv"),
    sep="\t",
    index=False
)
