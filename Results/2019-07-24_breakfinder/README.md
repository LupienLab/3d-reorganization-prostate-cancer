# Summary

This folder contains the breakpoints for structural variants called by `hic_breakfinder`.
`inter_expect_1Mb.hg38.txt` and `intra_expect_100kb.hg38.txt` were downloaded from the [Dixon Lab's Box archive](https://salkinstitute.app.box.com/s/m8oyv2ypf8o3kcdsybzcmrpg032xnrgx) on 2019-07-24 (at the time of writing, these files were last updated on 2018-05-22).

The breakpoint calls can be found in `Breakpoints/Min_1kbp/`.
`Breakpoints/Default/` contains breakpoints found by running `hic_breakfinder` without the `--min-1kb` option, since its runtime is much shorter.
