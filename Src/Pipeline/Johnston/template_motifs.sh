#!/bin/bash
# ===============================================================================================
# For Helix
#PBS -S /bin/bash
#PBS -N _job_type_-_sample_
#PBS -l walltime=8:00:00
#PBS -l mem=64gb
#PBS -l nodes=1:ppn=4
#PBS -q mainq

# Notes: mainq  - Max walltime = 7d. Max processors = 24. Mem = 256 Gb.
# Notes: fast	- Max walltime = 24h. Max processors = 16. Mem = 128 Gb.
# Notes: bigmem	- Max walltime = 24h. Max processors = 64. Mem = 2 Tb.
# Notes: gpu 	- Max walltime = 24h. Max processors = 16. Mem = 256 Gb. 2 NVIDIA Tesla K40m GPUs.
# ===============================================================================================

# ===============================================================================================
# For Synergy

#BSUB -J _job_type_-_sample_
#BSUB -n 4
#BSUB -R "span[hosts=1]"
#BSUB -W 2:00
#BSUB -o _job_type_-_sample__%J.out
#BSUB -e _job_type_-_sample__%J.err

# Notes: 56 cpu max
# ===============================================================================================

##########
# Setup
##########

# Need date for record-keeping and unique identifier
date="`date +%b%d_%T`"

# Take in values from make-submit + sed
topdir="_topdir_"
# topdir="/home/michael.johnston1/hiCAnalysis/hiC_523-567-583"
sample="_sample_"
# sample="G523_aggr"
job_type="_job_type_"
# job_type="arrowhead"
# outdir="_outdir_"
# # outdir=$topdir/postprocessing/arrowhead
logdir="_logdir_"
# logdir="$topdir/logs"
run_name="$job_type-$sample"

# Make this my working directory in case any loose files appear
cd $topdir

#############
# Directories for record-keeping
##############

# Record-keeping
# Redirect stdout ( > ) into a named pipe ( >() ) running "tee"
# See output on terminal
# As well as save it to file
exec > >(tee -i "$logdir"/$run_name.$date.log)
exec 2>&1

# ===========================================================

echo ""
echo topdir = $topdir
echo sample = $sample
echo job_type = $job_type
echo logdir = $logdir
echo outdir = $outdir
echo ""

##########
# Motifs
# Assign CTCF motifs to ends of Hiccups loop calls
# Ouputs: <loops>_with_motifs.txt
##########
# This one runs in ~1 min, but might be breaking the #CPU requested limit

echo ""
echo ""
echo "#===== Call Motifs ====="
echo Motifs Start: `date +%b%d_%T`

# Paths
hic_input="$topdir/juicer_out/$sample/aligned/inter_30.hic"

# Helix
# path_to_juicer_tools="/home/michael.johnston1/software/juicer_CPU/scripts/common/juicer_tools"
# motif_list="/home/michael.johnston1/reference/hg38.CTCF.motifs.txt"

# Synergy
path_to_juicer_tools="/home/mjohnston/software/juicer_CPU/scripts/common/juicer_tools"
# path_to_juicer_tools="/home/mjohnston/software/juicer_LSF/scripts/juicer_tools"
motif_list="/home/mjohnston/reference/hg38.CTCF.motifs.txt"

# ChIP data
# data_path="/home/michael.johnston1/hiCAnalysis/hiC_523-567-583/chip_data/gallo_ana/G523"
data_path=$topdir/chip_data/$sample

# Copy merged loops to allow data_name incorporation
# merged_loops="$TOPDIR/aligned/inter_30_loops.txt/merged_loops"
# merged_loops="$topdir/juicer_out/$sample/aligned/inter_30_loops_distantMod/merged_loops"
# merged_loops="$topdir/juicer_out/$sample/aligned/inter_30_loops_5-50-100/merged_loops"
merged_base="$topdir/juicer_out/$sample/aligned/inter_30_loops_distantMod/merged_loops"
merged_loops="${merged_base}.bedpe"
merged_loops_copy="${merged_base}_for_${sample}_ChIP.bedpe"
cp "$merged_loops" "$merged_loops_copy"

# Motifs
# motifs <genomeID> <bed_file_dir> <looplist> [custom_global_motif_list]
# File path to a directory which contains two folders: "unique" and "inferred".
# These folders should contain a combination of RAD21, SMC3, and CTCF BED files.
# If only CTCF data is available, use the same ChIP-Seq peaks in both the "unique" and "inferred" folders.

CMD="$path_to_juicer_tools motifs hg38 \
$data_path \
$merged_loops_copy \
$motif_list"

echo [CMD:] $CMD
$CMD

echo ""
echo "#===== Complete ====="
echo Motifs Done: `date +%b%d_%T`