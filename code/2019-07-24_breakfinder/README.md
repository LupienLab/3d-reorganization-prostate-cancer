# Summary

This folder contains the breakpoints for structural variants (SVs) called by `hic_breakfinder`.
This code can not be run on CodeOcean due to the reliance on patient Hi-C BAM files and privacy concerns.
Output from this folder, however, can be used on CodeOcean.

## Materials & Methods

[`hic_breakfinder`](https://github.com/dixonlab/hic_breakfinder/tree/30a0dcc6d01859797d7c263df7335fd2f52df7b8) was used to identify structural variants from the Hi-C data.

`inter_expect_1Mb.hg38.txt` and `intra_expect_100kb.hg38.txt` were downloaded from the [Dixon Lab's Box archive](https://salkinstitute.app.box.com/s/m8oyv2ypf8o3kcdsybzcmrpg032xnrgx) on 2019-07-24 (at the time of writing, these files were last updated on 2018-05-22).

## Results

### Structural variants can be detected with `hic_breakfinder`

`hic_breakfinder` identified between 6 and 78 SVs across all 12 samples, using default parameters.
The breakpoint calls can be found in [`../../results/2019-07-24_breakfinder/Breakpoints/Default/`](../../results/2019-07-24_breakfinder/Breakpoints/Default/).

Notably, the number of breakpoint pairs correspond to ERG over-expression.
The 7 ERG over-expressing samples are listed within the top 8 positions.
It is also notable that in only 1 sample (`PCa51852`) was the T2E fusion detected, despite the fusion being visible in the contact matrices of all 6 samples with that fusion.

### Differences in structural variant detection is not due to differences in sequencing depth

We hypothesized that the T2E fusion was not being detected in the remaining 5 samples due to differences in read depths, given that `PCa51852` happens to be the sample with the T2E fusion with the most filtered read pairs.

After downsampling `PCa51852` and re-running `hic_breakfinder` with default parameters, the T2E fusion was still detected, suggesting this was not an issue of sequencing depth.

### Differences in structural variant detection likely stems from detection thresholds

Our next hypothesis for why the T2E fusion was not detected in all 6 samples where it is known to exist was due to the detection threshold implicit in `hic_breakfinder`.
The closeness of the deletion to the diagonal of the matrix may not contrast enough with the local background in the first stage of submatrix finding, leading to it not being detected at the 100 kbp resolution.
Jesse Dixon created a branch of `hic_breakfinder` that allows a pre-set log-odds ratio threshold to be set with the `--thresh` option.
[`../../legacy/2019-07-24_breakfinder/Breakpoints/Low_Thresh_30/`](../../legacy/2019-07-24_breakfinder/Breakpoints/Low_Thresh_30/) contains breakpoints found by running `hic_breakfinder` with the `--thresh` option, to increase the sensitivity to breakpoint calls.

Fixing the threshold at 30, `PCa51852` is still the only sample where the T2E fusion is detected.
Fixing the threshold at 20, `PCa56413` has the T2E fusion detected.
Due to the size of the deleted region leading to the fusion, and the increased detection at a lower log-odds threshold, this suggests that the threshold for reasonably small SVs close to the diagonal is the limiting factor for detecting this SV.

### Refined structural variant calls do not identify new variants

To refine the detected SVs for accurate location, we re-ran `hic_breakfinder` with the `--min-1kb` option.

These detected breakpoints are all the same as the default parameters, with some of the SVs having their resolution reduced from 10 kbp to 1 kbp, except in 1 case in `PCa56413` where two SVs on chr14 were combined into a single event at 1 kbp resolution.

Overall, this refinement of SVs may provide a more accurate location for each detection, but does not identify new variants.

### Detected structural variants still need to be manually annotated

`hic_breakfinder` finds a submatrix that has a local maximum in signal, and returns the row and column coordinates for this submatrix.
It does not indicate what type of structural variant it has detected, however.
While Dixon _et al._, 2018, provide heuristics for determining what a given SV will appear like in the contact matrix (see Supplementary Figure 4), this is not incorporated into their statistical model.

To address this limitation, we manually curated each SV called by `hic_breakfinder` in each sample to resolve what type of SV each event was and its locus(i) to within 100 kbp.

| Structural variant type | Abbreviation[^1] |
| ----------------------- | ---------------- |
| Inversion               | `INV`            |
| Duplication             | `DUP`            |
| Insertion               | `INS`            |
| Deletion                | `DEL`            |
| Translocation           | `BND`            |

The manually curated SV calls can be found in [`../../results/2019-07-24_breakfinder/Breakpoints/Default/PCa*.breaks.manually-resolved.tsv`](../../results/2019-07-24_breakfinder/Breakpoints/Default/).

### A previously-observed T2E deletion is instead a translocation that affects chromatin topology at the insertion site

Despite the T2E fusion not being detected with default parameters in `PCa13848`, there was an event detected involving the T2E locus and chr14.

On the left is `PCa13848`, and on the right is another sample with the T2E fusion, but no evidence of a translocation (`PCa3023`).
chr21 is along the x-axis, and chr14 along the y-axis.
The insertion of this fragment from chr21 occurs near `chr14:35000000`.
The downstream end of the chr21 segment is inserted at the upstream end of the insertion site leading to _TMPRSS2_ being located next to _INSM2_ and _ERG_ being located next to _BRMS1L_.

Previously, this event was only detected as a deletion on chr21 using whole genome sequencing (see Supplementary Table S16 from Fraser et al., 2017.
Moreover, not only can this SV be reclassified from a deletion to a translocation, this event also affects the local chromatin topology at the insertion site.

The TAD in `PCa3023` at `chr14:35040000-35840000` is split into two TADs in `PCa13848` (`chr14:35040000-35720000` and `chr14:35720000-35840000`).

## Conclusions

We can robustly identify SV breakpoints from the Hi-C data using `hic_breakfinder`.
The breakpoints come in pairs and still need to be manually annotated, or combined into larger, complex events.
We see a difference in SV breakpoint pairs between T2E+ and T2E- patients.
Once complex SVs are constructed, we can return to this question.

## Footnotes

[^1]: These abbreviations are consistent with Delly, a commonly-used tool for detecting SVs from whole genome sequencing data.
