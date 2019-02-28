#!/bin/bash
# ===============================================================================================
#PBS -S /bin/bash
#PBS -N _job_type_-_sample_
#PBS -l walltime=48:00:00
#PBS -l mem=80gb
#PBS -l nodes=1:ppn=8

# Notes: fast	- Max walltime = 24h. Max processors = 16. Mem = 128 Gb.
# Notes: bigmem	- Max walltime = 24h. Max processors = 64. Mem = 2 Tb.
# Notes: gpu 	- Max walltime = 24h. Max processors = 16. Mem = 256 Gb. 2 NVIDIA Tesla K40m GPUs.
# ===============================================================================================

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

##########
# Arrowhead
# Contact Domains
# Ouputs: <size>_blocks
##########

echo ""
echo topdir = $topdir
echo sample = $sample
echo job_type = $job_type
echo logdir = $logdir
echo ""
echo "# ===== Call Arrowhead ====="
echo Time Start: `date +%b%d_%T`

# Path
path_to_juicer_tools="/home/michael.johnston1/software/juicer_CPU/scripts/common/juicer_tools"
hic_input="$topdir/juicer_out/$sample/aligned/inter_30.hic"
outdir="$topdir/juicer_out/$sample/aligned/inter_30_contact_domains.txt/"

# Make an output directory, if it doesn't exist
if [ -d  "$outdir" ]; then
   echo $outdir dir DOES exist
else
   echo $outdir dir DOES NOT exist -- Making dir
   mkdir "$outdir" || echo Unable to make diretory $outdir
fi

# Very limited test. Only Chr4. 25 kb res. Ignore sparse data
# juicebox arrowhead -c 4 -r 25000 --ignore_sparsity ../aligned/inter_30.hic ../postproc/arrowhead

echo \#=========100 000 blocks=========
$path_to_juicer_tools arrowhead -r 100000 $hic_input $outdir
echo \#=========50 000 blocks=========
$path_to_juicer_tools arrowhead -r 50000 $hic_input $outdir
echo \#=========25 000 blocks=========
$path_to_juicer_tools arrowhead -r 25000 $hic_input $outdir
echo \#=========10 000 blocks=========
$path_to_juicer_tools arrowhead -r 10000 $hic_input $outdir
echo \#=========5000 blocks===========
$path_to_juicer_tools arrowhead -r 5000 $hic_input $outdir

echo ""
echo "#===== Arrowhead Complete ====="
echo Arrowhead Done: `date +%b%d_%T`