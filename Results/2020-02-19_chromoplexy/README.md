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

## Results

### Complex SVs are common in prostate cancer and more frequent in _T2E_-fusion patients

Each connected component of the breakpoint graphs is a set of related breakpoints as a part of a complex event.
Below is a table with the number of connected components for each sample.

| Sample     | # Components | Largest number of connected breakpoints | T2E Status |
| ---------- | ------------ | --------------------------------------- | ---------- |
| `PCa13266` | 3            | 2                                       | No         |
| `PCa13848` | 17           | 6                                       | Yes        |
| `PCa14121` | 4            | 4                                       | No         |
| `PCa19121` | 7            | 4                                       | Yes        |
| `PCa3023`  | 15           | 12                                      | Yes        |
| `PCa33173` | 1            | 4                                       | No         |
| `PCa40507` | 16           | 10                                      | Yes        |
| `PCa51852` | 6            | 4                                       | Yes        |
| `PCa53687` | 7            | 12                                      | No         |
| `PCa56413` | 33           | 16                                      | Yes        |
| `PCa57054` | 15           | 12                                      | Yes        |
| `PCa57294` | 3            | 6                                       | No         |
| `PCa58215` | 4            | 6                                       | No         |

Complex events are components with > 2 connected breakpoints.
This occurs in 12/13 patients (only `PCa13266` does not have a complex event).
Like with the number of breakpoints detected, the T2E samples tend to have larger numbers of breakpoints involved in complex events.

## Conclusions

Merging redundant breakpoint calls demonstrates the amount of complex structural variants in these samples.
Moreover, it highlights the difference in mutational load between T2E+ and T2E- prostate cancer patients.
