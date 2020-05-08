# Summary

This folder contains the raw sequencing data from [Rhie _et al._, Nature Communications, 2019](https://doi.org/10.1038/s41467-019-12079-8).

The processed data available on [GEO Datasets](https://www.ncbi.mln.nig.gov/geo/query/acc.cgi?acc=GSE118629) is unfortunately aligned to the hg19 genome, so this pre-processing has to be redone to match with our data.

## Data

In the paper's methods, it describes that each samples was sequenced paired-end 75 bp to produce ~ 500 M read pairs per replicate using and Illumina HiSeq 2000.
The FASTQs coming from GEO do not show paired-end, but the sequences are 152 nts long.
So before aligning, these reads have to be split in the middle (producing two 76 bp sequences) and put into separate FASTQs.
These Hi-C samples were made with the MboI restriction enzyme, which is what we have used for our samples as well, so I don't need to create another artificially-digested genome with HiCUP.

## QC Checks

Looking at the FastQC reports for all samples, it appears that all samples look good, except for the high duplication rates in R2 for SRR8446386 and SRR8446385 (C42B replicates 1 and 2).
I'm unsure how to deal with these, given that the R1 mates for those samples have < 1% duplication rates.
So I'll have to treat these samples carefully, going forward.

