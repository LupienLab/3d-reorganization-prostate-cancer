# Summary

This folder contains results of using [TADsplimer](https://doi.org/10.1186/s13059-020-01992-7) on our Hi-C samples to compare TADs across samples, to assess whether structural variant breakpoints appear to split or merge TADs in a measurable way.

Due to TADsplimer's implementation, it cannot easily be installed as a Python/Conda package.
It's easiest point of entry is through a Docker container.
I've used Docker to perform these analyses, and the results should be reprodicuble, even if it doesn't fit the results/snakemake/compare format of the other `Results/` folders.
