# Summary

Upon earlier inspection of our data, we noticed that the _FOXA1_, _MYC_, and _AR_ genes were all located in close proximity to TAD boundaries in all 13 patients.
The _MYC_ locus is shown below.

![_MYC_ locus](../2019-10-24_higlass/Plots/MYC-locus.png)

These observations aligned with previous literature stating that housekeeping genes were also more frequently located near TAD boundaries [1] and that cohesin binds to DNA near actively transcribed genes [2].
These observations led us to the hypothesis that essential genes for the cell of interest are located near TAD boundaries (i.e. TADs form around these essential genes).

This folder contains the work that attempts to address that hypothesis.

## Materials & Methods

We use gene essentiality data from the [DepMap Project](https://depmap.org) for 5 prostate cancer cell lines (LNCaP, DU145, V16A, PC-3, and 22Rv1) to identify genes essential for prostate cancer.
We use the aggregated TAD boundaries from [`../2020-01-15_TAD-aggregation/`](../2020-01-15_TAD-aggregation/) from primary prostate cancer tumours.

## Results

### Genes are preferentially located near TAD boundaries

As seen below, the locations of TSSs are skewed towards TAD boundaries.
More than 50% of TSSs are located within the first 20% of a TAD, with a mode near 5%.

![Gene TSS locations across TADs](Plots/distance-density.png)
![Gene TSS locations across TADs CDF](Plots/distance-ecdf.png)

## Conclusions

## References

[1] Chen, Ke, Wu, Zhao, _et al._, Nature, 2019. doi: [10.1038/s41586-019-1812-0](https://doi.org/10.1038/s41586-019-1812-0)

[2] Busslinger _et al._, Nature, 2017. doi: [10.1038/nature22063](https://doi.org/10.1038/nature22063)