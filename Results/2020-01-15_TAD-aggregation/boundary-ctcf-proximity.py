"""
boundary-ctcf-proximity.py
==========

Calculate statistics related to how close LNCaP CTCF binding sites are to called TAD boundaries
"""

import os.path as path
import numpy as np
import pandas as pd
from tqdm import tqdm

# ==============================================================================
# Functions
# ==============================================================================
def distance(p: pd.Series) -> pd.Series:
    """
    Helper function for calculating the genomic distance between the TAD boundary and CTCF peak

    Parameters
    ==========
    p: pd.Series
        Row of the data table with `start_peak`, `start_bound`, and `end_peak` values
    """
    # return 0 if boundary is within the peak
    if (p.start_peak <= p.start_bound) and (p.start_bound <= p.end_peak):
        return 0
    # record whether the peak is upstream (-) or downstream (+) of the boundary
    # technically while the boundary could be between the start and end of the peak, the function will have already returned
    # 0, so we don't need to worry about that case here
    sgn = np.sign(p.start_peak - p.start_bound)
    return sgn * np.min([abs(p.start_peak - p.start_bound), abs(p.end_peak - p.start_bound)])



# ==============================================================================
# Data
# ==============================================================================
# load metadata
metadata = pd.read_csv(
    path.join("..", "..", "Data", "External", "LowC_Samples_Data_Available.tsv"),
    sep="\t",
    header=[0]
)
metadata = metadata.loc[metadata.Include == "Yes",:]
SAMPLES = ["PCa" + str(i) for i in metadata["Sample ID"]]

# load boundary-CTCF peak pairings for each sample
pairs = pd.concat(
    [
        pd.read_csv(
            path.join("CTCF", s + ".LNCaP-CTCF-peaks.bed"),
            sep="\t",
            header=None,
            names=[
                "chr_flank", "start_flank", "end_flank", "chr_bound", "start_bound", "end_bound", "persistence", "w",
                "chr_peak", "start_peak", "end_peak", "name_peak", "score_peak", "strand_peak", "signal", "p", "q", "peak"

            ],
            usecols=[
                "chr_bound", "start_bound", "end_bound", "persistence", "w", "chr_peak", "start_peak", "end_peak",

            ],
        ) for s in SAMPLES
    ],
    keys=SAMPLES,
    names=["SampleID"],
)
# convert SampleID index into a column
pairs.reset_index(inplace=True)
pairs.drop("level_1", axis=1, inplace=True)

# convert |-separated window size column into a proper list
pairs["w"] = pairs.w.str.split("|")

# create IDs for each boundary
pairs["Boundary_ID"] = pairs["chr_bound"].astype(str) + ":" + pairs["start_bound"].astype(str)

# ==============================================================================
# Analysis
# ==============================================================================
# calculate distance between CTCF peak and the TAD boundary
pairs["Distance"] = pairs.apply(distance, axis=1)

# create 5 kbp bins around the boundary
cuts = pd.interval_range(start=-200000 + 2500, end=200000 - 2500, freq=5000)
pairs["Distance_Bin"] = pd.cut(pairs["Distance"], cuts)

all_counts = pd.DataFrame({"Freq": pairs.groupby(["SampleID", "Distance_Bin"]).size() / pairs.groupby("SampleID").size()})
all_counts.reset_index(drop=False, inplace=True)

# ==============================================================================
# Save data
# ==============================================================================
# convert window size column into a comma-separated string
pairs["w"] = [",".join(pw) for pw in pairs["w"]]

# add locations of the bin for easier plotting
all_counts["Bin_Lower"] = [b.left for b in all_counts["Distance_Bin"]]
all_counts["Bin_Mid"] = [int((b.left + b.right) / 2) for b in all_counts["Distance_Bin"]]
all_counts["Bin_Upper"] = [b.right for b in all_counts["Distance_Bin"]]

# save to TSV
all_counts.to_csv(
    path.join("CTCF", "TAD-boundary.LNCaP-CTCF-peaks.distances.tsv"),
    sep="\t",
    header=True,
    index=False
)
