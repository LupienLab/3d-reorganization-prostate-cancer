# Summary

This folder contains work related to the altered histone acetylation surrounding structural variant breakpoints.

## Data

We use the H3K27ac ChIP-seq data, aligned to the hg38 reference genome, located in [`../../Data/Processed/2019-05-03_PCa-H3K27ac-peaks/`](../../Data/Processed/2019-05-03_PCa-H3K27ac-peaks/).
We also use the TADs identified in [`../2020-01-15_TAD-aggregation/`](../2020-01-15_TAD-aggregation/), and the test groups for these TADs listed in [`../2020-02-19_sv-disruption-TADs/sv-disruption-tests.TADs.tsv`](../2020-02-19_sv-disruption-TADs/sv-disruption-tests.TADs.tsv).

## Methods

To test for differential acetylation, we follow the procedure outlined in [DiffBind](https://www.bioconductor.org/packages/release/bioc/vignettes/DiffBind/inst/doc/DiffBind.pdf).
The analysis here, however, has a few features of structural variants and their relationship to TADs that must be considered:

1. The units under consideration are TADs, not peaks
2. None of the structural variant breakpoints are recurrent, aside from the _T2E_ fusion
3. Some overlapping TADs from different patients contain breakpoints resulting from different types of structural variants, and are considered separately

Like the RNA-seq analysis, these facts have the following implications for analysis:

1. Distinct features must be found that are consistent between all samples
2. All comparisons will have a "mutated" group of size $n = 1$ (except for the _T2E_ comparison)
3. Not all comparisons will have a "non-mutated" group of size $n = 11$; it depends on whether another patient has a breakpoint in an overlapping TAD

### Creating distinct features that are consistent across all samples

In canonical differential analysis for ChIP-seq data, peaks are called individually in each sample, and then peaks from all samples are merged to create a consensus peak set.
This is appropriate for ChIP-seq data since peaks are typically spaced apart, and peaks from the same sample are not typically merged together into the same consensus peak.
TADs, on the other hand, segment the entire genome and are not spaced far apart from each other.
This means that neighbouring TADs in the same sample will certainly be merged if the TAD boundaries in another sample are slightly different.
This means that large regions of the genome will almost certainly be grouped together, reducing the relevance of the final results, since each feature will span multiple neighbouring TADs.

To counter this behaviour, we take the induced segmentation approach described in [Zaborowski and Wilczyński, Journal of Computational Biology, 2019](https://doi.org/10.1089/cmb.2018.0162).
The genome is segmented into non-overlapping bins of variable width where the boundaries of the induced region correspond to any boundary in any of the samples being considered.

![Induced segmentation](Plots/induced-segmentation-schematic.png)

The resulting features are consistent across all samples, and because of the shared genomic coordinates, are also matched for other features like GC content and mappability.

### Normalization for non-ideal group sizes

Implication #2 is a typical non-starter for differential sequencing analysis; however, due to the nature of these structural variants, this is unavoidable.
To address this, we adopt a normalization technique described by [DESeq](https://bioconductor.org/packages/release/bioc/vignettes/DESeq/inst/doc/DESeq.pdf) (see "Working partially without replicates").
We use the "non-mutated" samples to calculate dispersion estimates for both groups using a `local` fit, while using scaling factors calculated from the mapped library size from each sample.
This relies on the assumption that most regions will not have differential acetylation between our "mutated" and "non-mutated" groups.
In our case, the smallest size of the "non-mutated" group is 8.
Note that this step is not taken for the _T2E_ fusion, since each group contains 6 samples and is appropriately balanced.

This approach is strictly forbidden in the `DiffBind` R package, so this analysis must be implemented in its own script, which we have done in [`altered-acetylation.R`](altered-acetylation.R).

### Combining results to return focus on TADs

After the differential analysis is performed, each induced segment has an array of associated statistics, including `baseMean`, `baseMean`, `log2FoldChange`, `lfcSE`, `stat`, `p`, and `padj`.
To relate the analysis back to the TAD from the mutated sample(s), the statistics are combined for all induced segments contained within the mutated sample(s)'s TAD.

![Reverse mapping of induced segments](Plots/reverse-map-induced-segments.png)

Statistics are combined as follows:

#### stat and p-value

p-values are combined using a weighted Stouffer's method, where the weights correspond to the fraction of the original TAD that the induced segment occupies.
Formally,

$$
z_i' = \frac{\sum_{j | o_j \cap D_i} w_j \cdot \phi^{-1}(1 - p_j)}{\sqrt{\sum_{j | o_j \cap D_i} w_j^2}}
$$

where $\phi$ is the standard normal CDF, $p_j$ is the p-value associated with induced segment $j$, and the weight for segment $j$ is

$$
w_j = \frac{|o_j \cap D_i|}{|D_i|}
$$

The resulting p-value for TAD $i$, is provided by a two-sided z-test: $p_i' = 2 \cdot \phi(-|z_i'|)$

#### baseMean, log2FoldChange, and lfcSE



#### Multiple testing corrections

We are only testing for differential acetylation in the TADs containing structural variants, this is not a genome-wide comparison.
Thus, we only correct for the number of tests performed per normalization group.
