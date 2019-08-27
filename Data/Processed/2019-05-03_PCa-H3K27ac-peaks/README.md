# Summary

This folder contains called peaks from all the CPC-GENE samples processed by Ken and Stanley.
Since these FASTQs originate from a collection of flowcells, they're aggregated into this single folder for uniform processing.
See the folders for specific flowcells to see how they were aligned.

The ENCODE Blacklist [1] was downloaded from https://www.encodeproject.org/annotations/ENCSR636HFF/.

> Note: Since MACS requires Python 2.7 and all other analyses in this project use Python 3, a separate conda environment needed to be created to accommodate this tool.
> You can find the YAML file for this environment in `macs.yaml`.
> Use this environment before running the peak calling in this folder.
> In the `Snakefile`, the `conda` keyword is used, and the `--use-conda` flag is used in the command-line execution of this folder.

## Results



## References

1. Amemiya, H. M., Kundaje, A. & Boyle, A. P. The ENCODE Blacklist: Identification of Problematic Regions of the Genome. Scientific Reports 9, 9354 (2019).

