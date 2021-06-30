# Summary

This folder is part of the reviewer response.
It was suggested that we use an alternative method for broadly comparing tumour and benign samples with the help of tools like [HiCRep](https://github.com/TaoYang-dev/hicrep).
Here we address this comment.

## Methods

Instead of the original HiCRep R package, we use the Python implementation of the algorithm from the [hicreppy](https://pypi.org/hicreppy) package.
We do this because our Hi-C files are stored in the Cooler format, and this new package natively works with these files.
We also do this because when we initially tried using the R package, we ran into memory constraints from the algorithm itself.
These constraints appear to be addressed in the Python package.

