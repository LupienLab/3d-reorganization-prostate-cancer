# Summary

This folder contains the hierarchical domains called by TopDom for all 13 prostate samples.
We us TopDom [1] (in `topdom.R`) to call topologically associated domains (TADs) using a set of window size parameters to give a hierarchical view of domain organization.

## Results

Using a contact matrix at a resolution of 40 kbp, we use a set of values for the window parameter ($w \in [3, 30]$) to detect TADs across a variety of length scales.

We then aggregate these calls for each window size into a hierarchy TAD calls based on their size and location.

## References

[1]: Shin, H. et al. TopDom: an efficient and deterministic method for identifying topological domains in genomes. Nucleic Acids Res 44, e70–e70 (2016). doi: [10.1093/nar/gkv1505](https://doi.org/10.1093/nar/gkv1505)
