#!/bin/bash
# ===============================================================================================
#PBS -S /bin/bash
#PBS -N _job_type_-_sample_
#PBS -l walltime=2:00:00
#PBS -l mem=40gb
#PBS -l nodes=1:ppn=1
#PBS -q mainq

# Notes: mainq  - Max walltime = 168h. Max processors = 24. Mem = 256 Gb.
# Notes: fast	- Max walltime = 24h.  Max processors = 16. Mem = 128 Gb.
# Notes: bigmem	- Max walltime = 24h.  Max processors = 64. Mem = 2 Tb.
# Notes: gpu 	- Max walltime = 24h.  Max processors = 16. Mem = 256 Gb. 2 NVIDIA Tesla K40m GPUs. Only one job at a time.
# ===============================================================================================

##########
# Setup
##########

# # # Load modules
# module load cuda

# =========================================================

##########
# Setup
##########

# # Load cuda module
# # Only for hiccups
# module load cuda/7.5.18

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
# Calculate Resolution
# Find how small bind can be made while still maintaining minimal coverage
# Ouputs: aligned/eigen_<chr>_<bp>.txt
##########

echo ""
echo ""
echo "#===== Calculate Resolution ====="
echo Time Start: `date +%b%d_%T`

# Usage: calculate_map_resolution.sh <merged_nodups file> <50bp coverage file>
#   <merged_nodups file>: file created by Juicer containing all valid+unique read pairs
#   <50bp coverage file>: where to write the 50bp coverage file; if this file is non-empty, the 50bp coverage vector wont be recalculated
# Took ~ 40 minutes per sample

# # From Rao 2014:
# # The ‘‘map resolution’’ is the smallest locus size such that 80%
# # of loci have at least 1,000 contacts.

script_to_run="/home/michael.johnston1/software/juicer_CPU/scripts/calculate_map_resolution.sh"
input_merge_nodups="$topdir/juicer_out/$sample/aligned/merged_nodups.txt"
bp_cov_file="$topdir/juicer_out/$sample/aligned/50bp_coverage_vector.txt"

# Check if output exists
# Move if it does
if [ -s  $bp_cov_file ]
	then
		echo Coverage file: $bp_cov_file DOES exist.
		# mv "$bp_cov_file" "$bp_cov_file".bak || echo Unable to move
		echo ***Warning: coverage file already exisits***
		echo May need to move / delete if you want to recalculate.
		echo ""
	else
		echo Coverage file: $bp_cov_file DOES NOT exist -- Proceeding
		echo ""
fi

# Run
CMD="$script_to_run $input_merge_nodups $bp_cov_file > $topdir/juicer_out/$sample/aligned/resolution.txt"
echo ""
echo ""
echo [CMD:] $CMD
eval $CMD

echo ""
echo "Run end time: `date +%b%d_%T`"