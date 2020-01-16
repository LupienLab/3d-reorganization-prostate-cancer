# Summary

In [`../2019-07-08_TADs/`](../2019-07-08_TADs/), we identify TADs at a variety of length scales by using different smoothing window sizes, $w$.
In this folder, we present a method for aggregating these TAD calls across length scales to resolve boundaries that are nearby the same and to establish the "order" of a TAD boundary (i.e. how many sub-TADs have this locus as a boundary).

## Materials & Methods

Starting with TAD calls from `TopDom` over a set of window sizes, $W = \{w\}$, we identify all boundaries for all TADs.

### Boundary aggregation

Boundaries are then aggregated together over $W$, keeping track of which window size each boundary was identified with.
This leaves a set of boundaries, where each boundary corresponds to a set of window sizes, $V \subset W$.
We define the "order" of a particular boundary as $o(b) = |V|$.

### Resolving conflicting boundaries

Due to the stochastic nature of sequencing methods and TAD calling as performed by `TopDom`, some boundaries from different window sizes will be adjacent to each other, but fundamentally represent the same boundary.
To resolve these adjacent boundary calls, we identify "conflicting boundaries" as boundaries that are:

1. within 1 bin (e.g. 40 kbp) of each other
1. called at different window sizes

Any set of conflicting boundaries are resolved to a single boundary as follows.
Let $B = \{b_i\}$ be a set of conflicting boundaries, and let $V_i$ be the set of window sizes corresponding to the boundary $b_i$.

### When $|B| = 2$

If $|B| = 2$ and all of the window sizes for boundary $b_i$ are greater than all of those for $b_{j \ne i}$ (i.e. $u > v \forall u \in V_i, v \in V_j$), then $b_i$ is selected as the "resolved boundary", $\hat{b}$, $b_j$ is discarded, and we set the windows sizes of $\hat{b}$ to be $\hat{V} = V_i \cup V_j$.

If, however, $u \ngeq v \forall u \in V_i, v \in V_j$, then the "resolved boundary" is set to whichever boundary has the highest order (i.e. $\hat{b} = \argmax_{b \in B} o(b)$, $\hat{V} = \cup_i V_i$).

### When $|B| = 3$

The symmetry of this scenario suggests that taking the mean of these boundaries is likely accurate, on average, so we select the "resolved boundary" to be the middle of the three (i.e. $\hat{b} = median(B)$, $\hat{V} = \cup_i V_i$)).

### When $|B| > 3$

In this scenario, multiple bins are called as a TAD boundary across window sizes, without a unique resolved boundary.
To computational simplicity, like in the case of $|B| = 2$ without a clear ordering of window sizes, we designate the "resolved boundary" to be the boundary with the largest order (i.e. $\hat{b} = \argmax_{b \in B} o(b)$, $\hat{V} = \cup_i V_i$).

### TAD reconstruction

With the resolved boundaries across all window sizes, we then reconstruct TADs at each window size by segmenting the genome according to the location of each resolved boundary.
All of the above steps are written as code in [`aggregate-TADs.R`](aggregate-TADs.R).

## Results

