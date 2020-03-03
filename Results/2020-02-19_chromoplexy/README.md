# Summary

When observing the contact matrices, we noted how many events seemed linked together by their location.
This phenomenon of multiple SVs where chromosomal ends are shuffled without leading to copy number losses is termed "chromoplexy", and has previously been observed as particularly common in prostate cancer \Cref{Li2020,Baca2013}.
Here, we explore this phenomenon using our Hi-C data.

## Data

We use the called breakpoints from `hic_breakfinder` in [`../2019-07-24_breakfinder/Breakpoints/Default/`](../2019-07-24_breakfinder/Breakpoints/Default/).
We also make use of the RNA-seq data from these same samples \Cref{Chen2019}, and the TADs we previously called in [`../2020-01-15_TAD-aggregation/`](../2020-01-15_TAD-aggregation/).

## Methods

### Graph construction of SV breakpoints

Similar in design to the ChainFinder algorithm previously described \Cref{Baca2013}, we represent each breakpoint as a node in a graph.
Each row of the `hic_breakfinder` output contains a pair of breakpoints corresponding to the bounding coordinates of the aberrant submatrix.
These two nodes are connected via an edge, coloured according to the breakpoint type assigned to it by manual annotation (see the `Type` column in [`../2019-07-24_breakfinder/Breakpoints/Default/PCa*.breaks.sorted.manually-resolved.tsv`](../2019-07-24_breakfinder/Breakpoints/Default/)).
Breakpoints are subsequently connected with an edge if they are within 100 kbp of each other.
This tolerance distance was chosen due to the granularity of the breakpoint calls, since each breakpoint pair is identified at a contact matrix resolution of 10 kbp or 100 kbp (most breakpoints called at 1 Mbp resolution appear to be a false positive due to effects of compartmentalization).
This produces a graph of SV breakpoints for each patient, where every connected component of the graph (i.e. sets of inter-connected breakpoints) is a complex event.

### Detecting SV breakpoints that alter local chromatin topology

We consider each breakpoint end from each SV detected.
For each breakpoint at locus $[s, e)$ in a given patient, we look for other patients that have a breakpoint end within 500 kbp (i.e. overlapping $[s - \delta, e + \delta), \delta = 500 000$) to identify similarly mutated patients.
A local version of the BPscore \Cref{Zaborowski2019} is calculated for the identified TADs within this $[s - \delta, e + \delta)$ window for each pair of samples.
Using these calculations, we use two approaches for determining whether there is an alteration of the local topology; one binary classification method, and one permutation test method.

#### Binary classification of local topology

Patients are then assigned into one of two groups using hierarchical clustering with the matrix of pairwise BPscore values as a distance matrix.
If the clustering perfectly stratifies the mutated samples from the non-mutated samples (i.e. the clustering matches the mutation status in this locus), then the local topology is counted as altered as a result of the SV.

#### Permutation test for differences in local topology



### Hypothesis testing for differences in RNA abundance

Conventional methods for differential gene expression, such as DESeq2 \Cref{Love2014}, EdgeR \Cref{Robinson2010}, and Sleuth \Cref{Yi2018} require replicates for each condition being tested.
For our case of comparing a sample with an SV to samples without, recurrent events are rare, thus leading to 1-vs-many comparisons, which is insufficient for these previous methods.

To address this shortcoming, we developed a different null hypothesis testing framework by aggregating genes in TADs containing the SV breakpoints (implemented in [`test-dge.py`](test-dge.py)).

1. For each breakpoint, identify the overlapping TAD(s).
  In the case of an insertion, it may be a single TAD; in the case of a large deletion, it may be many.
2. Identify all genes lying within these affected TADs, according to the GENCODE v33 reference \Cref{GENCODE}.
3. Perform a z-transformation on the abundance values (in FPKM) for these genes, given by
  $$
  z_i = \frac{x_{i,mut} - \mu_{i,non-mut}}  {\sigma_{i,non-mut}}
  $$
  Under the null hypothesis, that these linked breakpoints do not affect the expression of the nearby genes within the same TADs, the random variable for the expression of gene is distributed like the non-mutant samples, namely
  $$
  \mathbb{E}[X_{i,mut}] = \mu_{i,non-mut}
  $$
  $$
  \text{Var}[X_{i,mut}] = \sigma^2_{i,non-mut}
  $$
  Thus, each gene's z-score should come from a distribution with mean 0 and standard deviation 1.

3. For each set of related genes, perform a two-sided t-test on these z-scores.

## Results

### Chromplexy is common and more frequent in _T2E_-fusion patients

### Breakpoints rarely alter the local chromatin topology

Using the method described above for detecting changes to local topology as a result of and SV breakpoint, we find that a small minority of SVs alter the local topology (11/462 breakpoint ends, 2.38%).

SVs that do alter the local topology include:

* a translocation of the deleted _TMPRSS2_-_ERG_ locus being inserted into chr14 in `PCa13848`
* multiple complex events on chr4 of `PCa3023`
* a translocation of a deleted segment on chr12 and inserted into chr17 in `PCa3023` near _NCOR1_, _TTC19_, and _SNORD163_
* one end of a tandem duplication on chr3 in `PCa53687` by _GAP43_
* a complex, indeterminate event on chr15 in `PCa53687`
* an apparent chromosome arm swap between chr7 and chr19 in `PCa53687`
* multiple chained events on chr3 of `PCa56413`
* a duplication on chr10 of `PCa56413`

Only one breakpoint of a detected SV appeared to alter the local topology.
There was no event detected where both breakpoints in a pair altered the local topology.
If SVs are altering gene expression, it is likely not through establishing or altering TAD boundaries, but by interfering with _cis_-regulatory interactions through other means.

### Chromoplexy alters the expression of genes within TADs containing breakpoints

Using the method described above, we identified $n$ chromoplexic events that significantly altered gene expression of genes within the same TADs as the linked breakpoints.

![Distribution of expression fold changes for each complex event](Plots/sv-disruption.fold-change.png)

While most genes did not have a large change to their expression, we observed that approximately one third (33.1%) of all genes had fold changes greater than 2.

![Chromoplexy leads to altered expression in 1/3 of nearby genes](Plots/sv-disruption.fold.ecdf.png)

To reduce false positives stemming from lowly-expressed genes (thus producing large fold-changes and z-scores, despite not having large absolute changes in expression), we consider the absolute change in RNA abundance.
We see that the majority of genes with large fold changes (i.e. $|\log_2(\text{fold change})| \ge 1$) have small absolute differences in read abundances between the mutated and non-mutated samples (i.e. $|\bar{x}_{mut} - \bar{x}_{non-mut}|$ < 0.1 FPKM).

![Fold change vs absolute RNA abundance difference between mutated and non-mutated samples](Plots/sv-disruption.fold-change-vs-difference.png)

This reduces the number of differentially expressed genes to 605 (from the original 3882 without setting a threshold on absolute RNA abundance difference).

## Conclusions

