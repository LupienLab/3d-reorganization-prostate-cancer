# Summary

This folder is a revamped version of the ideas explored in [`../2020-10-23_t2e-loops/`](../2020-10-23_t2e-loops/).
Here, our goal is to test whether gene expression changes observed in T2E+/- groups can be explained by epigenetic changes in H3K27ac at enhancers and their looping contact with gene promoters.
Our scientific hypothesis is that these changes can be explained as such, but that given the high false negative rate for loop calls, as well as the size of the loops, that we will only be able to provide sufficient evidence for a small number of differentially expressed genes.
The null hypothesis is that there is no relationship between the differential expression of genes and these loops, which we will test.

## Data

We will use a few different datasets:

1. GENCODE gene annotations for promoter regions
2. Catalogue of H3K27ac peaks across the 12 primary tumour samples (the entire list is contained in [`../2020-06-12_sv-disruption-acetylation/Acetylation/T2E/t2e.all.tsv`](../2020-06-12_sv-disruption-acetylation/Acetylation/T2E/t2e.all.tsv))
3. TADs identified in [`../2020-08-29_TADs-downsampled/Aggregated-TADs/separated-TADs/`](../2020-08-29_TADs-downsampled/Aggregated-TADs/separated-TADs/) (window size = 20)
4. Catalogue of loops for all 17 primary samples, as well as the T2E-specific and non-T2E-specific loops identified in [`../2020-10-23_t2e-loops/loops.T2E-specific.tsv`](../2020-10-23_t2e-loops/loops.T2E-specific.tsv) and [`../2020-10-23_t2e-loops/loops.nonT2E-specific.tsv`](../2020-10-23_t2e-loops/loops.nonT2E-specific.tsv), respectively.

## Methods

### GRN construction

We first construct a gene regulatory network (GRN) for each annotated gene.
We do this by identifying the parent TAD for each gene, then creating a graph where the gene promoter and each enhancer within that TAD is a node.
An edge connects two elements that overlap each other's genomic coordinates (such as the peak at a gene's promoter and the promoter region itself) or if there is a loop connecting the two anchors together.

### GRN filtering for T2E-specific epigenetic changes

We will test for an association between differential expression and the presence of differential loops and acetylation changes at loop anchors.
We consider all genes whose GRNs match the following criteria:

1. $\exists \ge 1$ differential loop in the GRN between T2E+ and T2E-
2. The differential loops are all gained or all lost in the T2E+ group
3. $\exists \ge 1$ differential peak in the loop anchors that are all in the same direction as each other, and in the same direction as the differential loop
4. GRN does not have an opposing peak change elsewhere in the GRN

Filters 3 and 4 are to reduce the ambiguous effects of multiple simultaneous epigenetic changes, and to limit the effect of false differential loop classification by prioritizing H3K27ac ChIP-seq data, which has stronger support for it.

The construction of these GRNs can be found in `grn.py`.

### Hypothesis testing

After filtering down to a set of testable genes, we are left with 3 groups of genes:

1. Genes with a gain of activating epipgenetic changes in T2E+
2. Genes with a loss of activating epigenetic changes in T2E+ (equivalently, a gain of activiating epigenetic change in T2E-)
3. Genes with ambiguous or no epigenetic changes to their GRN

We test against the null hypothesis that the three groups have the same mean fold change in expression between TE+ and T2E- samples.
This test requires an ANOVA to test for a difference in means between these 3 groups.

## Results
