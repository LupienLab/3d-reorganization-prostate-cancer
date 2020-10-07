# Summary

In this folder, we consider the loop calls in primary prostate samples, both benign and tumour.

## Data

We use the multi-resolution cooler files from each Hi-C experiment stored in [`../../Data/Processed/2019-06-18_PCa-LowC-sequencing/`](../../Data/Processed/2019-06-18_PCa-LowC-sequencing/).

## Methods

Loops were previously called using cLoops.
We also call loops using `cooltools call-dots`, which considers the enrichment of contacts over an expected background.
Results of the loop calls can be found in [`../../Data/Processed/2019-06-18_PCa-LowC-sequencing/Loops/`](../../Data/Processed/2019-06-18_PCa-LowC-sequencing/Loops/).

## Results

`cooltools call-dots` only return a few hundred loops per sample, which seems cery incorrect, to me.
There should be thousands of loops, given that the current working model of gene regulation involves gene promoters looping to enhancer regions.

After switching to [Mustache](https://doi.org/10.1186/s13059-020-02167-0), we identified thousands of loops per sample, which appears much more reasonable.
