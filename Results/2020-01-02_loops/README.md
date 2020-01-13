# Summary

This folder contains data related to calling loops from the Low-C data using [cLoops](https://github.com/YaqiangCao/cLoops).

## Results

### Called loops

cLoops requires setting 2 parameters (as it is based on DBSCAN): _eps_ (the distance within which two paired-end tags are classified as neighbours) and _minPts_ (the minimum required number of paired-end tags required to be called a loop).
Running the default command with cLoops allows the program to select a minimum self-ligation distance cutoff after using a set of _eps_ and _minPts_ parameters for the candidate loop detection.
The distance cutoff for each sample, DBSCAN parameters, and the number of loops called are listed in `Loops/loop-stats.tsv` (`Distance_Cutoff` in units of bp).

![Loops called across all samples](Plots/loop-stats.png)

### Loops overlapping with active cis-regulatory elements

We hypothesize that these loops correspond to promoter-enhancer loops and TAD boundaries, while loops that do not correspond to these features are false detections by `cLoops`.
We start by seeing how many loops have 0, 1, or 2 endpoints overlapping active cis-regulatory elements (CREs; as defined by H3K27ac peaks).
We hypothesize that:

* loops involving 2 active CREs are promoter-enhancer loops (where both are active in prostate cancer)
* loops involving 1 active CREs are promoter-enhancer loops (where the enhancer is poised)
* loops involving 0 active CREs are TAD boundaries

![Loop anchors overlapping active cis-regulatory elements](Plots/loop-CRE-overlap.sum.proportion.png)


This suggests that many of the loops being called are in fact loops involved in gene regulation.