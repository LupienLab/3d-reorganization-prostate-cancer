#!/bin/bash
# ===============================================================================================
#PBS -S /bin/bash
#PBS -N juicer-SAMP-N_early
#PBS -l walltime=24:00:00
#PBS -l mem=2000gb
#PBS -l nodes=1:ppn=63
#PBS -q bigmem

# Notes: fast	- Max walltime = 24h. Max processors = 16. Mem = 128 Gb.
# Notes: bigmem	- Max walltime = 24h. Max processors = 64. Mem = 2 Tb.
# Notes: gpu 	- Max walltime = 24h. Max processors = 16. Mem = 256 Gb. 2 NVIDIA Tesla K40m GPUs.
# ===============================================================================================

##########
# Setup
##########

DATE="`date +%b%d_%T`"
echo "Run start time: `date +%b%d_%T`"

#############
# Directories for record-keeping
##############

# Assume this is being run from the top directory
# TOPDIR=`pwd`
TOPDIR=/gpfs/home/mjohnston/data_gpfs/hiC_523-567-583_redo/juicer_out/G523_aggr

# Make an output directory, if it doesn't exist
RUN_NAME="juicer-SAMP-N"
DIR_NAME="out_logs"
OUT_DIR="$TOPDIR/$DIR_NAME"

if [ -d  "$OUT_DIR" ]; then
   echo $DIR_NAME dir DOES exist
else
   echo $DIR_NAME dir DOES NOT exist -- Making dir
   mkdir "$OUT_DIR" || echo Unable to make diretory $DIR_NAME
fi

# Record-keeping
# Redirect stdout ( > ) into a named pipe ( >() ) running "tee"
# See output on terminal
# As well as save it to file
exec > >(tee -i "$OUT_DIR"/$RUN_NAME.$DATE.log)
exec 2>&1

# ===========================================================

##########
# Juicer call begins here
##########
# Directories of required inputs and references
##########




echo ""
echo ""
echo ===== Call Juicer =====
echo Time Start: `date +%b%d_%T`

# # Define topdir
# TOPDIR="/home/michael.johnston1/hiCAnalysis/many_fastq"

# Location of Juicer manager script
JUICER_SCRIPT="/export/home/mjohnston/software/juicer_LSF/scripts/juicer.sh"

# JUICEDIR is a directory that contains Juicer elements
# Including scripts, references, restriction enzyme sites and more
# JUICEDIR="/home/michael.johnston1/hiCAnalysis/juice_Dir"
JUICEDIR="/export/home/mjohnston/software/juicer_LSF"

# Number of threads to use
THREADS=16

# Genome
GENOME_ID="hg38"

# Restriction Site
RES_SITE='MboI'

# About / Description
ABOUT_DESC='SAMP-N'

# Stage
# [stage]: must be one of "merge", "dedup", "final", "postproc", or "early".
#     -Use "merge" when alignment has finished but the merged_sort file has not
#      yet been created.
#     -Use "dedup" when the files have been merged into merged_sort but
#      merged_nodups has not yet been created.
#     -Use "final" when the reads have been deduped into merged_nodups but the
#      final stats and hic files have not yet been created.
#     -Use "postproc" when the hic files have been created and only
#      postprocessing feature annotation remains to be completed.
#     -Use "early" for an early exit, before the final creation of the stats and
#      hic files

# Add in below
# -S $STAGE \

# STAGE='early'
# STAGE='merge'
# STAGE='postproc'


# # Remove pesky existing folders during troubleshooting
# # DO NOT
# # BE DUMB
# # AND OVERWRITE
# # YOUR WORK
# if [ -d "$TOPDIR/aligned" ]; then
# 	rm -r $TOPDIR/aligned
# fi

# All jobs to normal queue on Synergy
short_queue="normal"
long_queue="normal"

# Run juicer
# Make sure -z and -p options come before -g

# Usage: juicer.sh [-g genomeID] [-d topDir] [-q queue] [-l long queue] [-s site]
#                  [-a about] [-R end] [-S stage] [-p chrom.sizes path]
#                  [-y restriction site file] [-z reference genome file]
#                  [-C chunk size] [-D Juicer scripts directory]
#                  [-Q queue time limit] [-L long queue time limit] [-r] [-h] [-x]

CMD="$JUICER_SCRIPT \
-d $TOPDIR \
-z $JUICEDIR/references/hg38.fasta \
-p $JUICEDIR/references/hg38.chrom.sizes \
-q $short_queue \
-l $long_queue \
-g $GENOME_ID \
-s $RES_SITE \
-D $JUICEDIR \
-a $ABOUT_DESC \
-y $JUICEDIR/restriction_sites/hg38_MboI.txt"

echo ""
echo [CMD:] $CMD
echo ""
$CMD

echo ""
echo "Run end time: `date +%b%d_%T`"