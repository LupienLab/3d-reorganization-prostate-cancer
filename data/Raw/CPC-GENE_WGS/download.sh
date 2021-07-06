#!/bin/sh
#SBATCH -p download
#SBATCH -t 7-00:00:00
#SBATCH -c 4
#SBATCH --job-name CPCG_down

sftp -P 13579 lupienlab@cpcgene.fraserlab.ca:/feec27a8c8141e70fa960f0df5552a7e.CPCG_0336_Pr_P_PE_572_WG_130927_h239_0210_BC2DA7ACXX_4_NoIndex_R1.fastq.gz FASTQs/
lupienlab
