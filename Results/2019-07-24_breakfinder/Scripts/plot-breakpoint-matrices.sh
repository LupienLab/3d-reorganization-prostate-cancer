#!/bin/bash
# ==================================================================================================
# Environment
# ==================================================================================================
PLOTTING_SCRIPT="../2019-07-08_TADs/plot-tads.py"
COOL_DIR="../../Data/Processed/2019-06-18_PCa-LowC-sequencing/Contacts"
WINDOW=2500000

# declarative array (aka hash table) for chromosome sizes
declare -A CHROM_SIZES=( \
    ["chr1"]=248956422 \
    ["chr2"]=242193529 \
    ["chr3"]=198295559 \
    ["chr4"]=190214555 \
    ["chr5"]=181538259 \
    ["chr6"]=170805979 \
    ["chr7"]=159345973 \
    ["chr8"]=145138636 \
    ["chr9"]=138394717 \
    ["chr10"]=133797422 \
    ["chr11"]=135086622 \
    ["chr12"]=133275309 \
    ["chr13"]=114364328 \
    ["chr14"]=107043718 \
    ["chr15"]=101991189 \
    ["chr16"]=90338345 \
    ["chr17"]=83257441 \
    ["chr18"]=80373285 \
    ["chr19"]=58617616 \
    ["chr20"]=64444167 \
    ["chr21"]=46709983 \
    ["chr22"]=50818468 \
    ["chrX"]=156040895 \
    ["chrY"]=57227415 \
    ["chrM"]=16569 \
)

# ==================================================================================================
# Command line interface
# ==================================================================================================
usage () {
    echo "Create contact matrix plots including highlighted breakpoints."
    echo "Usage: ./$(basename $0) -i INPUT [-o OUTDIR] [-h]"
    echo ""
    echo "Options:"
    echo "-h, --help                Show help"
    echo "-i, --input=INPUT         Breakpoints called by hic_breakfinder"
    echo "-o, --outdir=OUTDIR       Output directory for images"
}

while test $# -gt 0; do
    case "$1" in
        -h|--help)
            usage
            exit 0
            ;;
        -i)
            shift
            if test $# -gt 0; then
                export DATA=$1
            else
                echo "No input data specified"
                exit 1
            fi
            shift
            ;;
        --input*)
            export DATA=`echo $1 | sed -e 's/^[^= ]*[= ]//g'`
            shift
            ;;
        -o)
            shift
            if test $# -gt 0; then
                export OUTDIR=$1
            else
                echo "No outdir specified"
                exit 1
            fi
            shift
            ;;
        --outdir*)
            export OUTDIR=`echo $1 | sed -e 's/^[^=]*=//g'`
            shift
            ;;
        *)
            usage
            exit 1
            ;;
    esac
done

if [ $OUTDIR == "" ]; then
    export OUTDIR="."
fi

# ==================================================================================================
# Execution
# ==================================================================================================
# iterate through each breakpoint in the BEDPE file
n_breakpoints=$(wc -l ${DATA} | cut -f 1 -d " ")
i=0

while IFS=$'\t' read -r -a row; do
    i=$(($i + 1))
    echo "${i} of ${n_breakpoints}" >&2
    # parse line from data file
    chrom1="${row[0]}"
    start1="${row[1]}"
    end1="${row[2]}"
    chrom2="${row[3]}"
    start2="${row[4]}"
    end2="${row[5]}"
    sample="${row[6]}"

    # plotting regions to see breakpoint in its context
    #   1 Mbp on either side, if possible
    col_lo=$(( $start1 < $WINDOW ? 0 : $start1 - $WINDOW ))
    col_hi=$(( $end1 + $WINDOW > ${CHROM_SIZES[$chrom1]} ? ${CHROM_SIZES[$chrom1]} : $end1 + $WINDOW ))
    row_lo=$(( $start2 < $WINDOW ? 0 : $start2 - $WINDOW ))
    row_hi=$(( $end2 + $WINDOW > ${CHROM_SIZES[$chrom2]} ? ${CHROM_SIZES[$chrom2]} : $end2 + $WINDOW ))

    # format strings for higlight regions
    hc="${chrom1}:${start1}-${end1}"
    hr="${chrom2}:${start2}-${end2}"

    # format strings for display region
    reg1="${chrom1}:${col_lo}-${col_hi}"
    reg2="${chrom2}:${row_lo}-${row_hi}"

    # create output file
    output="${OUTDIR}/${sample}.${chrom1}_${start1}_${end1}.${chrom2}_${start2}_${end2}.png"
    python ${PLOTTING_SCRIPT} --zmin=-4.5 --zmax=-1.5 -o ${output} --hc ${hc} --hr ${hr} -r2 ${reg2} ${COOL_DIR}/${sample}.mcool::/resolutions/40000 ${reg1}

    # check that script completed successfully
    echo -e "${sample}.${chrom1}_${start1}_${end1}.${chrom2}_${start2}_${end2}\t$?"
done < ${DATA}
