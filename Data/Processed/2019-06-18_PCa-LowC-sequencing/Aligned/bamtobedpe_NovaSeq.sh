#!/bin/bash
#SBATCH -t 2-00:00:00
#SBATCH --mem=6000
#SBATCH -J BamtoBedPE
#SBATCH -p himem
#SBATCH --array=0-12
#SBATCH -c 1
#SBATCH -N 1
#SBATCH -o ./log/%x-%j.out


echo 'Hello!'

## ID for task
idx=$SLURM_ARRAY_TASK_ID

echo ${idx}

input_files=(PCa40507.agg.hicup.name-sorted.bam PCa51852.agg.hicup.name-sorted.bam PCa13266.agg.hicup.name-sorted.bam PCa53687.agg.hicup.name-sorted.bam PCa13848.agg.hicup.name-sorted.bam PCa56413.agg.hicup.name-sorted.bam PCa57054.agg.hicup.name-sorted.bam PCa14121.agg.hicup.name-sorted.bam PCa57294.agg.hicup.name-sorted.bam PCa19121.agg.hicup.name-sorted.bam PCa58215.agg.hicup.name-sorted.bam PCa3023.agg.hicup.name-sorted.bam PCa33173.agg.hicup.name-sorted.bam)

input=${input_files[${idx}]}

output_files=(PCa40507.agg.hicup.name-sorted.bedpe PCa51852.agg.hicup.name-sorted.bedpe PCa13266.agg.hicup.name-sorted.bedpe PCa53687.agg.hicup.name-sorted.bedpe PCa13848.agg.hicup.name-sorted.bedpe PCa56413.agg.hicup.name-sorted.bedpe PCa57054.agg.hicup.name-sorted.bedpe PCa14121.agg.hicup.name-sorted.bedpe PCa57294.agg.hicup.name-sorted.bedpe PCa19121.agg.hicup.name-sorted.bedpe PCa58215.agg.hicup.name-sorted.bedpe PCa3023.agg.hicup.name-sorted.bedpe PCa33173.agg.hicup.name-sorted.bedpe) 

output=${output_files[${idx}]}

echo "Sample name: "${input}

module load bedtools

bedtools bamtobed -i ${input} -bedpe > ${output}

echo 'End!'
