# ==============================================================================
# Meta
# ==============================================================================
# plot-t2e-insertion.py
# ------------------------------------------------
# Description: Plot neo-loops and neo-TADs that result from the T2E insertion in PCa31848
# Author: James Hawley

import argparse
PARSER = argparse.ArgumentParser(
	formatter_class=argparse.ArgumentDefaultsHelpFormatter
)
PARSER.add_argument(
	"sample",
	type=str,
	help="Sample ID to plot",
	choices=[
		"PCa13266",
		"PCa13848",
		"PCa14121",
		"PCa19121",
		"PCa3023",
		"PCa33173",
		"PCa40507",
		"PCa51852",
		"PCa53687",
		"PCa56413",
		"PCa57294",
		"PCa58215",
	]
)
ARGS = PARSER.parse_args()

import os.path as path
import pandas as pd
import cooler
import matplotlib.pyplot as plt
import mpl_toolkits
# from mpl_toolkits.axes_grid1 import colorbar
from neoloop.visualize.core import *


# ==============================================================================
# Data
# ==============================================================================
DIR = {
	"res": path.join("..", "..", "results", "2021-07-26_neoloop"),
	"H3K27ac": path.join("..", "..", "data", "Processed", "2019-05-03_PCa-H3K27ac-peaks", "Tracks"),
}
DIR["cooler"] = path.join(DIR["res"], "cooler")
DIR["loop"] = path.join(DIR["res"], "neoloops")
DIR["tad"] = path.join(DIR["res"], "neotads")
DIR["assembly"] = path.join(DIR["res"], "assembled-SVs")
DIR["plot"] = path.join(DIR["res"], "Plots")

# load cooler file
clr = cooler.Cooler(path.join(DIR["cooler"], ARGS.sample + ".res_10000bp.cnv-corrected.cool"))

# load list of oncogenes from the NeoLoopFinder package
oncogene_list = [line.rstrip() for line in open("oncogenes.txt")]

# load rearranged assembly
assemblies_df = pd.read_csv(
	path.join(DIR["assembly"], ARGS.sample + ".res_10000bp.assemblies.txt"),
	sep="\t",
	header=None,
	names=["Assembly", "Junction", "Breakpoint_1", "Breakpoint_2"],
)


# ==============================================================================
# Analysis
# ==============================================================================
# process each assembly in the table into something well-formatted that
# neoloop.visualize.core.Triangle will parse
assemblies = {
	a[1]: "\t".join([a[1], a[2], a[3], a[4]]) for a in assemblies_df.itertuples()
}


# ==============================================================================
# Plots
# ==============================================================================
# plot each assembly
for assembly_name, assembly_str in assemblies.items():
	vis = Triangle(
		clr,
		assembly_str,
		n_rows=5,
		figsize=(7, 5.2),
		track_partition=[5, 0.8, 0.8, 0.2, 0.5],
		correct='weight',
		span=300000,
		space=0.08,
	)
	vis.matrix_plot(vmin=0, cbr_fontsize=9)
	vis.plot_chromosome_bounds(linewidth=2)
	# plot RNA-seq signal
	# vis.plot_signal('RNA-Seq', 'enc_SCABER_RNASeq_rep1.bw', label_size=10, data_range_size=9, max_value=0.5, color='#E31A1C')
	# plot H3K27ac signal
	vis.plot_signal(
		"H3K27ac",
		path.join(DIR["H3K27ac"], ARGS.sample + "_FE.sorted.filtered.bw"),
		label_size=10,
		data_range_size=9,
		max_value=20,
		color='#6A3D9A'
	)
	# plot genes
	vis.plot_genes(filter_=oncogene_list, fontsize=10)
	# plot chromosome names
	vis.plot_chromosome_bar(name_size=13, coord_size=10)
	# add loops and neo-loops
	vis.plot_loops(
		path.join(DIR["loop"], ARGS.sample + ".res_10000bp.neoloops.tsv"),
		face_color="none",
		marker_size=40,
		cluster=True,
		onlyneo=False,
	)
	# save the figure
	vis.outfig(
		path.join(DIR["plot"], ".".join([ARGS.sample, "res_10000bp", assembly_name, "png"])),
		dpi=300
	)
	vis.outfig(
		path.join(DIR["plot"], ".".join([ARGS.sample, "res_10000bp", assembly_name, "pdf"])),
		dpi=300
	)
