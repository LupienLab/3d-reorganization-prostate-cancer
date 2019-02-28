#!/bin/bash
# ===============================================================================================
#PBS -S /bin/bash
#PBS -N _job_type_-_sample_
#PBS -l walltime=8:00:00
#PBS -l mem=128gb
#PBS -l nodes=1:ppn=16
#PBS -q gpu

# Notes: fast	- Max walltime = 24h. Max processors = 16. Mem = 128 Gb.
# Notes: bigmem	- Max walltime = 24h. Max processors = 64. Mem = 2 Tb.
# Notes: gpu 	- Max walltime = 24h. Max processors = 16. Mem = 256 Gb. 2 NVIDIA Tesla K40m GPUs.
# ===============================================================================================

##########
# Setup
##########

# Load cuda module
# Only for hiccups
module load cuda/7.5.18

# =========================================================

# Need date for record-keeping and unique identifier
date="`date +%b%d_%T`"
# Take in values from make-submit + sed
topdir="_topdir_"
# topdir="/home/michael.johnston1/hiCAnalysis/hiC_523-567-583"
sample="_sample_"
# sample="G523_aggr"
job_type="_job_type_"
# job_type="arrowhead"
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
echo ""

##########
# Hiccups
# Loops
# Ouputs: inter_30_loops.txt/merged_loops
##########
# HiCCUPS
# juicebox hiccups [-m matrixSize] [-k normalization (NONE/VC/VC_SQRT/KR)] 
# [-c chromosome(s)] [-r resolution(s)] [-f fdr] [-p peak width] [-i window] 
# [-t thresholds] [-d centroid distances] [--ignore_sparsity]
# <hicFile> <outputDirectory> [specified_loop_list]

# Use --ignore_sparsity if low resolution maps

echo "# ===== Call HiCCUPS ====="
echo Time Start: `date +%b%d_%T`

# Path
path_to_juicer_tools="/home/michael.johnston1/software/juicer_CPU/scripts/common/juicer_tools"
hic_input="$topdir/juicer_out/$sample/aligned/inter_30.hic"

# # Modification for GM12878
# echo "Using link for GM12878"
# hic_input="https://hicfiles.s3.amazonaws.com/hiseq/gm12878/combined_30.hic"

# Default parameters
#######
# echo Running with defaults
# outdir="$topdir/juicer_out/$sample/aligned/inter_30_loops"
# CMD="$path_to_juicer_tools hiccups -m 2048 -r 5000,10000,25000 $hic_input $outdir"

# Distant Loop Modification
#######
# Parameters used here are
# HICCUPS defaults for 5,10,25 kb
# Taken from Rao Cohesin Loss paper for 50 and 100 kb
echo Running with Distant Loop Modification
outdir="$topdir/juicer_out/$sample/aligned/inter_30_loops_distantMod"
CMD="$path_to_juicer_tools hiccups -m 2048 \
-r 5000,10000,25000,50000,100000 \
-f 0.1,0.1,0.1,0.1,0.1 \
-p 4,2,1,2,1 \
-i 7,5,3,4,2 \
-d 20000,20000,50000,100000,200000 \
$hic_input \
$outdir"

echo [CMD:] $CMD
$CMD

echo ""
echo "#===== HiCCUPS Complete ====="
echo HiCCUPS Done: `date +%b%d_%T`