# Summary

This folder contains digested genomes for various restriction enzymes.

## Characteristics of fragments for various enzymes

All of the following are for hg19.

### HindIII

| Statistic | Value        |
| --------- | ------------ |
| count     | 8.376010e+05 |
| mean      | 3.594164e+03 |
| std       | 4.868967e+04 |
| min       | 6.000000e+00 |
| 25%       | 9.190000e+02 |
| 50%       | 2.275000e+03 |
| 75%       | 4.671000e+03 |
| max       | 3.021669e+07 |

### MboI

| Statistic | Value        |
| --------- | ------------ |
| count     | 7.127583e+06 |
| mean      | 4.223964e+02 |
| std       | 1.657816e+04 |
| min       | 4.000000e+00 |
| 25%       | 1.260000e+02 |
| 50%       | 2.650000e+02 |
| 75%       | 5.430000e+02 |
| max       | 3.000220e+07 |

## For updating information

To add a new restriction enzyme:

* Open `config.tsv`
* Add a new line with the restriction motif and a caret (^) indicating the cut site

> Note: these are used implicitly by HiCUP.
> The start sites in the `.txt` files are 1-indexed, not 0.

## Conclusions

MboI is a much more frequent cutter than HindIII, with more binding sites and smaller fragment sizes (both on average and with a smaller variance).
Therefore, it's a better cutter to go with, and may help with reducing duplication rates.
