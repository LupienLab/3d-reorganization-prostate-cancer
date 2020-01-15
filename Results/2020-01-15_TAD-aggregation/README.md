# Summary

In [`../2019-07-08_TADs/`](../2019-07-08_TADs/), we identify TADs at a variety of length scales by using different smoothing window sizes, $w$.
In this folder, we present a method for aggregating these TAD calls across length scales to resolve boundaries that are nearby the same and to establish the "order" of a TAD boundary (i.e. how many sub-TADs have this locus as a boundary).

## Materials & Methods

Starting with TAD calls from `TopDom` over a set of window sizes, $W = \{w\}$, 