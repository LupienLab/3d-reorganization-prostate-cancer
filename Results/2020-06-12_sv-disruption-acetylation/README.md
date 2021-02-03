# Summary

This folder contains work related to the altered histone acetylation surrounding structural variant breakpoints.

## Materials

We use the H3K27ac ChIP-seq data, aligned to the hg38 reference genome.
We also use the TADs previously identified, and the test groups for these TADs listed in the compiled SVs folder.

## Methods

To test for differential acetylation of H3K27ac consensus peaks in individual samples identified by their structural variants., we follow the procedure outlined in [DiffBind](https://www.bioconductor.org/packages/release/bioc/vignettes/DiffBind/inst/doc/DiffBind.pdf).

### Normalization for non-ideal group sizes

The lack of recurrence in SVs between patients is unavoidable and leads to unbalanced experimental designs.
To address this, we adopt a normalization technique described by [DESeq](https://bioconductor.org/packages/release/bioc/vignettes/DESeq/inst/doc/DESeq.pdf) (see "Working partially without replicates").
We use the "non-mutated" samples to calculate dispersion estimates for both groups using a `local` fit, while using scaling factors calculated from the mapped library size from each sample.
This relies on the assumption that most regions will not have differential acetylation between our "mutated" and "non-mutated" groups.
In our case, the smallest size of the "non-mutated" group is 8.
Note that this step is not taken for the _T2E_ fusion, since each group contains 6 samples and is appropriately balanced.

This approach is strictly forbidden in the `DiffBind` R package, so this analysis must be implemented in its own script, which we have done in [`altered-acetylation.R`](altered-acetylation.R).

#### Multiple testing corrections

We are only testing for differential acetylation in the TADs containing structural variants, this is not a genome-wide comparison.
Thus, we only correct for the number of tests performed per normalization group.
