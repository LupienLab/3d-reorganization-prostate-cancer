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
from genomic_interval import GenomicInterval, overlapping, find_tad

# ==============================================================================
# Constants
# ==============================================================================
# promoter region offset (upstream, downstream)
PROM_OFFSET = (1500, 500)

# get chromosome sizes
hg38 = nc.get_chrominfo("hg38")
CHROM_SIZES = hg38.chrom_lengths

# ==============================================================================
# Functions
# ==============================================================================

# ==============================================================================
# Data
# ==============================================================================
# load GENCODE gene annotations
genes = pd.read_csv(
    path.join(
        "..", "..", "Data", "External", "GENCODE", "gencode.v33.all-genes.bed"
    ),
    sep="\t",
    header=None,
    names=["chr", "start", "end", "strand", "gene_id", "gene_name"],
)

# load catalogue of H3K27ac peaks and differential test results
peaks = pd.read_csv(
    path.join(
        "..",
        "2020-06-12_sv-disruption-acetylation",
        "Acetylation",
        "T2E",
        "t2e.all.tsv",
    ),
    sep="\t",
    header=[0],
)

# load catalogue of loop calls from all 17 samples
loops = pd.read_csv(
    path.join("..", "2020-10-06_loops", "Loops", "merged-loops.sample-counts.tsv"),
    sep="\t",
    header=[0],
)

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

# ==============================================================================
# Analysis
# ==============================================================================
# create placeholder for the GRN
G = {gid: nx.Graph() for gid in genes["gene_id"]}

# add the promoter and H3K27ac peaks within the TAD to the graph for each gene
for gene in genes.itertuples():
    prom = GenomicInterval(
        gene.chr,
        gene.strand == "+" ? max(0, gene.start - PROM_OFFSET[0]) : min(CHROM_SIZES[gene.chr], gene.start + PROM_OFFSET),
        gene.end,
        {
            "row": gene.Index,
            "id": gene.gene_id,
            "name": gene.gene_name,
            "type": "promoter",
        },
    )
    # add the promoter to the relevant GRN
    G[gene.gene_id].add_node(prom)
    # get parent TAD for this gene

# ==============================================================================
# Plots
# ==============================================================================

# ==============================================================================
# Main
# ==============================================================================


def main():
    """
    Main
    """
    pass


if __name__ == "__main__":
    PARSER = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    ARGS = PARSER.parse_args()
    main()
