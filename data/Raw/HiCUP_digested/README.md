# Summary

This folder contains digested genomes for various restriction enzymes.

## Characteristics of fragments for various enzymes

All of the following are for hg38.

### HindIII

| Statistic | Value        |
| --------- | ------------ |
| count     | 8.619840e+05 |
| mean      | 3.513807e+03 |
| std       | 4.290735e+04 |
| min       | 6.000000e+00 |
| 25%       | 8.880000e+02 |
| 50%       | 2.216000e+03 |
| 75%       | 4.603000e+03 |
| max       | 3.021669e+07 |

### MboI/DpnII

| Statistic | Value        |
| --------- | ------------ |
| count     | 7.203658e+06 |
| mean      | 4.206982e+02 |
| std       | 1.463577e+04 |
| min       | 4.000000e+00 |
| 25%       | 1.270000e+02 |
| 50%       | 2.670000e+02 |
| 75%       | 5.480000e+02 |
| max       | 3.000290e+07 |


## For updating information

To add a new restriction enzyme:

* Open `config.tsv`
* Add a new line with the restriction motif and a caret (^) indicating the cut site

> Note: these are used implicitly by HiCUP.
> The start sites in the `.txt` files are 1-indexed, not 0.

## Conclusions

MboI is a much more frequent cutter than HindIII, with almost 10 times more binding sites and 10 times smaller fragment sizes (both on average and with a smaller variance).
Therefore, it's a better cutter to go with, and may help with reducing duplication rates.
