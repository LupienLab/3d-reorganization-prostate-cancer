# Summary

Upon earlier inspection of our data, we noticed that the _FOXA1_, _MYC_, and _AR_ genes were all located in close proximity to TAD boundaries in all 13 patients.
The _FOXA1_ locus is shown below.

![_FOXA1_ locus](../2019-10-24_higlass/Plots/FOXA1-locus.png)

These observations aligned with previous literature stating that housekeeping genes were also more frequently located near TAD boundaries.
These observations led us to the hypothesis that essential genes for the cell of interest are located near TAD boundaries (i.e. TADs form around these essential genes).

This folder contains the work that attempts to address that hypothesis.

## Materials & Methods

We use gene essentiality data from the [DepMap Project](https://depmap.org) for 5 prostate cancer cell lines (LNCaP, DU145, V16A, PC-3, and 22Rv1) to identify genes essential for prostate cancer.
We use the aggregated TAD boundaries from [`../2020-01-15_TAD-aggregation/`](../2020-01-15_TAD-aggregation/) from primary prostate cancer tumours.

## Results

### Genes are uniformly distribution across TADs

![qq-plot of gene TSS location across TADs](Plots/PCa3023.proximity.qqplot.png)

## Conclusions
