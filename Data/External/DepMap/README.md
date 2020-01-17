# Summary

This folder contains gene essentiality data from the [DepMap Project](https://depmap.prg), which measures the essentiality of genes in different cell lines using siRNA and CRISPR knockout screens.
We use the 19Q4 data for both the Achilles dataset (CRISPR knockout) [1-3] and RNA interference [1, 4].

## DepMap Data Description

According to the DepMap portal (as of 2020-01-17), these are 10 cell lines derived from prostate cancer cells.
To keep only these data, we filter the genetic dependency scores with `filter-depmap.R`.
There is only 1 prostate cancer cell line in the CRISPR-Cas9 data (VCaP) and 7 in the RNAi data (22Rv1, DU145, LNCaP, MDA PCa 2b, NCI-H660, PC3, and VCaP).

For ease of comparison, we keep our analysis to the RNAi data.

# References

[1] DepMap, Broad (2019): "DepMap 19Q4 Public". figshare. Dataset [doi:10.6084/m9.figshare.9201770.v2](https://doi.org/10.6084/m9.figshare.9201770.v2).

[2] Robin M. Meyers, Jordan G. Bryan, James M. McFarland, Barbara A. Weir, ... David E. Root, William C. Hahn, Aviad Tsherniak. "Computational correction of copy number effect improves specificity of CRISPR-Cas9 essentiality screens in cancer cells". _Nature Genetics_ (2017). [doi:10.1038/ng.3984](https://doi.org/10.1038/ng.3984).

[3] Dempster, J. M., Rossen, J., Kazachkova, M., Pan, J., Kugener, G., Root, D. E., & Tsherniak, A. "Extracting Biological Insights from the Project Achilles Genome-Scale CRISPR Screens in Cancer Cell Lines". _bioRxiv_, [doi:10.1101/720243](https://doi.org/10.1101/720243).

[4] James M. McFarland, Zandra V. Ho, Guillaume Kugener, Joshua M. Dempster, Phillip G. Montgomery, Jordan G. Bryan, John M. Krill-Burger, Thomas M. Green, Francisca Vazquez, Jesse S. Boehm, Todd R. Golub, William C. Hahn, David E. Root, Aviad Tsherniak. "Improved estimation of cancer dependencies from large-scale RNAi screens using model-based normalization and data integration". _Nature Communications_ (2018). [doi:10.1038/s41467-018-06916-5](https://doi.org/10.1038/s41467-018-06916-5)
