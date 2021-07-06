#!/bin/bash
#SBATCH -p all
#SBATCH --mem=12G
#SBATCH -t 5-00:00:00
#SBATCH --job-name encrypt_tumour_hic_batch2

module load java

java -jar ../../../../ega-cryptor-2.0.0.jar -i FASTQs/ -o FASTQs/encrypted

