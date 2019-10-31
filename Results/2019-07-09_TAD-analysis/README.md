# Summary

This folder contains the analyses and interpretation of the TAD calls on all 13 prostate cancer samples done in `../2019-07-08_TADs/`.

## Results

### TAD similarity

To see whether TADs are called similarly in different samples, we calculate which TADs are consistent in pairwise comparisons.
We do this by counting TADs called in separate samples with at least 60% bi-directional overlap (i.e. ${(a, b) : a,b \in A \cap B, \frac{|a \cap b|}{|b|} \ge 0.6, \frac{|a \cap b|}{|a|} \ge 0.6}$, where $A$ and $B$ are the regions in both samples).

We find that on average, ~68% of TADs are consistently called in other samples for small TADs ($w = 3$), whereas ~ 86% of TADs are consistently called for larger TADs (> 500 kbp, $w \ge 10$).

![TAD consistency](Plots/tad-similarity-counts.png)

### Finding patients subsets with differential TADs

Our first approach to determining sets of TADs that are unique to a particular subset of patients.
