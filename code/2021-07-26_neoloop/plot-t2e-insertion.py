# ==============================================================================
# Meta
# ==============================================================================
# plot-t2e-insertion.py
# ------------------------------------------------
# Description: Plot neo-loops and neo-TADs that result from the T2E insertion in PCa31848
# Author: James Hawley

import os.path as path
import cooler
from neoloop.visualize.core import *

# ==============================================================================
# Data
# ==============================================================================
DIR = {
	"res": path.join("..", "..", "results", "2021-07-26_neoloop"),
}
DIR["cooler"] = path.join(DIR["res"], "cooler")
DIR["assembly"] = path.join(DIR["res"], "assembled-SVs")
DIR["plot"] = path.join(DIR["res"], "Plots")

# load cooler file
clr = cooler.Cooler(path.join(DIR["cooler"], "PCa13848.res_10000bp.cnv-corrected.cool"))

# load list of oncogenes from the NeoLoopFinder package
oncogene_list = [line.rstrip() for line in open("oncogenes.txt")]

# load rearranged assembly
assemblies = [line.rstrip() for line in open(path.join(DIR["assembly"], "PCa13848.assemblies.txt")]

# ==============================================================================
# Plots
# ==============================================================================
# assembly of interest is the C2 assembly on chr21
vis = Triangle(
	clr,
	assembly,
	n_rows=5,
	figsize=(7, 5.2),
	track_partition=[5, 0.8, 0.8, 0.2, 0.5],
	correct='weight',
	span=300000,
	space=0.08
)
vis.matrix_plot(vmin=0, cbr_fontsize=9)
vis.plot_chromosome_bounds(linewidth=2)
# plot RNA-seq signal
vis.plot_signal('RNA-Seq', 'enc_SCABER_RNASeq_rep1.bw', label_size=10, data_range_size=9, max_value=0.5, color='#E31A1C')
# plot H3K27ac signal
vis.plot_signal('H3K27ac', 'SCABER_H3K27ac_pool.bw', label_size=10, data_range_size=9, max_value=20, color='#6A3D9A')
# plot genes
vis.plot_genes(release=75, filter_=List, fontsize=10)
# plot chromosome names
vis.plot_chromosome_bar(name_size=13, coord_size=10)
# save the figure
vis.outfig(
	path.join(DIR["plot"], "PCa13848.t2e-insertion.png"),
	dpi=300
)

