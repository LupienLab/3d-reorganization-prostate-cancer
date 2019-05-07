# Summary

Using the H3K27ac peaks (re-aligned to hg38) for these primary samples, we'd like to know what resolution is required to distinguish peaks and often have them in separate bins, on average.

## Results

### H3K27ac peaks characterization

![Peak counts](Stats/peak-counts.png)

![Peak sizes](Stats/peak-sizes.png)

Note the different orders in the two previous plots.

The distance between peaks, which is needed to understand how high of a resolution we'll need for our contact matrices if we want to explore the interactions between putative CREs.

![Distance between peaks](Stats/peak-dists.png)

The median distance ranges between 238 and 501 bp, with means skewed much higher.

![Distance between peaks](Stats/peak-dists-centres.png)

Sticking with the mean distance, we would need to get to 5 kbp resolution to have 1 putative CRE per bin.
This should be accomplished with > 300 M filtered read pairs, as seen from Table S2 in [Rao _et al._](https://doi.org/10/xqj).

The medians show that this resolution may not be ideal, since the distances between peaks is skewed up by large distances.

## Conclusions

A resolution of 5 kbp is below the mean distance between H3K27ac peaks, but getting to a resolution below 1 kbp is desired to have a high likelihood that each putative CRE would belong to a single bin.
