# Summary

This folder uses and pre-processes raw Hi-C data in multiple benign prostate cell lines.
These come from 2 main papers: PrEC (3 replicates) from [Taberlay _et al._, Genome Research, 2014](https://doi.org/10.1101/gr.201517.115) and RWPE1 (2 replicates) from [Rhie _et al._, Nature Communications, 2019](https://doi.org/10.1038/s41467-019-12079-8).
Details are listed in [`config.tsv`](config.tsv).

As listed in the methods for each paper, the RWPE1 samples were digested using the MboI restriction enzyme and sequenced with paired-ends at 76 bp.
The PrEC samples were digested using and sequenced with BglII paired-ends at 101 bp.
The FASTQs downloaded from GEO don't come in pairs, so they need to be split before pre-processing.

