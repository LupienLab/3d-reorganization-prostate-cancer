#!/bin/bash
#SBATCH -p all
#SBATCH --mem=12G
#SBATCH -t 5-00:00:00
#SBATCH --job-name encrypt_benign_hic_batch1

module load java

java -jar ../../../../ega-cryptor-2.0.0.jar -i FASTQs/ -o FASTQs/encrypted

