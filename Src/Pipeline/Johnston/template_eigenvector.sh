#!/bin/bash
# ===============================================================================================
#PBS -S /bin/bash
#PBS -N eigen-_SAMP-N_-_chr-name_
#PBS -l walltime=24:00:00
#PBS -l mem=20gb
#PBS -l nodes=1:ppn=1
#PBS -q bigmem

# Notes: fast	- Max walltime = 24h. Max processors = 16. Mem = 128 Gb.
# Notes: bigmem	- Max walltime = 24h. Max processors = 64. Mem = 2 Tb.
# Notes: gpu 	- Max walltime = 24h. Max processors = 16. Mem = 256 Gb. 2 NVIDIA Tesla K40m GPUs.
# ===============================================================================================

# Note that for this post-processing, only the HiCCUPs step requires GPUs
# However, the large TMP storage is also required
# The GPU node can only handle one job at a time
# So may as well request all CPUs available

##########
# Setup
##########

# =========================================================

########
# Make TMP drive directories on bigmem/gpu
########

# Need date for record-keeping and unique identifier
DATE="`date +%b%d_%T`"
echo "Run start time: `date +%b%d_%T`"

# This is the 20 Tb temporary storage disk on bigmem/gpu
TMP="/local_scratch"

# Receive directory names from make-submit script
TOPDIR="$TMP/_directory-names_"
# TOPDIR="$TMP/michael.johnston1_Jul31_16:44:53_Sample_A_G523_P13"

SAMPLE=`basename $TOPDIR | cut -d "_" -f 3-`
# SAMPLE=`basename $TOPDIR | cut -d "_" -f 5-`

# Make directories if necessary
DIR_NAME=$TOPDIR
if [ -d  $DIR_NAME ]; then
   echo $DIR_NAME dir DOES exist
else
   echo $DIR_NAME dir DOES NOT exist -- Making dir
   mkdir $DIR_NAME || echo Unable to make diretory $DIR_NAME
fi

# Make this my working directory in case any loose files appear
cd $TOPDIR
echo ""
echo TOPDIR $TOPDIR
echo ""


#############
# Directories for record-keeping
##############

# Make an output log directory, if it doesn't exist
RUN_NAME="eigen-_SAMP-N_-_chr-name_"
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

##########
# Eigenvector
# Sign changes indicate switching between A / B compartments
# Ouputs: aligned/eigen_<chr>_<bp>.txt
##########

echo ""
echo ""
echo '#===== Call Eigenvector =====#'
echo Eigenvector Start: `date +%b%d_%T`

# Path
path_to_juicer_tools="/home/michael.johnston1/software/juicer_CPU/scripts/common/juicer_tools"

# Make directories if necessary
DIR_NAME="aligned/eigenvector"
OUT_DIR="$TOPDIR/$DIR_NAME"

if [ -d  "$OUT_DIR" ]; then
   echo $DIR_NAME dir DOES exist
else
   echo $DIR_NAME dir DOES NOT exist -- Making dir
   mkdir "$OUT_DIR" || echo Unable to make diretory $DIR_NAME
fi

# Eigenvector
# eigenvector <NONE/VC/VC_SQRT/KR> <hicFile(s)> <chr> <BP/FRAG> <binsize> [outfile]
# Force finer resolution with "p"
# Note that the sign does not indicate which compartment
# Sign flips indicate compartment switching 
# Compartment will need to be determined for each chromosome

# Do substitutions using sed from make-submit script
CMD="$path_to_juicer_tools eigenvector KR -p\
      $TOPDIR/aligned/inter_30.hic \
      _chr-name_ \
      BP _resolution_ \
      $OUT_DIR/eigen__chr-name___resolution_bp.txt"
      echo [CMD:] $CMD
      $CMD

# Unnecessary intermediate step now that eignevalueflipping works

# echo '#=== Converting to wig ===#'
# sed \
# "1i track type=wiggle_0 name=eigen__chr-name___resolution_bp description=eigen-sign color=255,0,0 altColor=0,255,0 autoScale=on\
# \nfixedStep chrom=chr_chr-name_ start=1 step=_resolution_ span=_resolution_"\
# $OUT_DIR/eigen__chr-name___resolution_bp.txt \
# > $OUT_DIR/eigen__chr-name___resolution_bp.wig

echo ""
echo Eigenvector Done: `date +%b%d_%T`