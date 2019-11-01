# Summary

This folder contains the analyses and interpretation of the TAD calls on all 13 prostate cancer samples done in `../2019-07-08_TADs/`.

## Results

### TAD similarity

To see whether TADs are called similarly in different samples, we calculate which TADs are consistent in pairwise comparisons.
We do this by considering TADs called in a given sample and counting whether that region is called as a TAD in at least one other sample with at least 60% bi-directional overlap (i.e. ${(a, b) : a,b \in A \cap B, \frac{|a \cap b|}{|b|} \ge 0.6, \frac{|a \cap b|}{|a|} \ge 0.6}$, where $A$ and $B$ are the regions in both samples).

We find that on average, ~68% of TADs are consistently called in other samples for small TADs ($w = 3$), whereas ~ 86% of TADs are consistently called for larger TADs (> 500 kbp, $w \ge 10$).
The general trend can be seen below.

![TAD consistency](Plots/tad-similarity-counts.png)

We see that for $w \ge 10$, the percentage of TADs in a given sample that are called similarly in another sample remains close to 85% on average, and is consistent across window sizes.
There is a monotonic reduction in mean similar TAD percentage for $w < 10$ for all samples.
These results suggest that large scale organization (on lengths scales of 400 kbp or larger) are largely similar across all samples, whereas small scale organization (< 400 kbp) differentiates samples, and may be where individual differences in etiology occur.

### Finding patient subsets with differential TADs

Our first approach to determining sets of TADs that are unique to a particular subset of patients
