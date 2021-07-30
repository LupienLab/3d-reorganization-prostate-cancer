# Summary

This folder contains code for a secondary attempt to identify novel enhancer-hijacking and neo-loops being formed around structural variants.

## Data

We use the contact matrices from all 17 primary samples in cooler format.

## Methods

We make use of the [NeoLoopFinder](https://github.com/XiaoTaoWang/NeoLoopFinder) scripts to identify multiple features of each sample.

### Identifying CNVs

We first identify putative CNVs with the `calculate-cnv` script.
We need to run the first instance of `calculate-cnv` on an internet-accessible cluster partition, since it downloads a GC content bigwig file from Northwestern University's Box account.
If we use the same cache for all subsequent `calculate-cnv` calls, then it doesn't need to download this file again.
According to the Supplementary Information of NeoLoopFinder's paper, it should take ~ 1h to compute CNVs for a single Cooler file, using < 2G of memory.

### Using structural variant breakpoints to create assemblies


### Calling neo-loops and neo-TADs

### Plotting

