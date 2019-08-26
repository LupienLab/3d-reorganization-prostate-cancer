# Summary

This folder contains the breakpoints for structural variants called by `hic_breakfinder`.
`inter_expect_1Mb.hg38.txt` and `intra_expect_100kb.hg38.txt` were downloaded from the [Dixon Lab's Box archive](https://salkinstitute.app.box.com/s/m8oyv2ypf8o3kcdsybzcmrpg032xnrgx) on 2019-07-24 (at the time of writing, these files were last updated on 2018-05-22).

The breakpoint calls can be found in `Breakpoints/Min_1kbp/`.
`Breakpoints/Default/` contains breakpoints found by running `hic_breakfinder` without the `--min-1kb` option, since its runtime is much shorter.
`Breakpoints/Low_Thresh/` contains breakpoints found by running `hic_breakfinder` without the `--thresh 30` option, to increase the sensitivity to breakpoint calls.
This is done in an attempt to find the T2E fusion.

## Results

### Lower threshold for detection

Despite setting the threshold lower at 30, PCa51852 remains the only sample where the T2E fusion is explicitly called.
Again, it is visible in all samples in the contact matrix, but is not completely detectable from this algorithm, we suspect due to its proximity.
The fusion occurs only within 3 Mbp, so it is difficult to distinguish these from the surrounding contacts.
