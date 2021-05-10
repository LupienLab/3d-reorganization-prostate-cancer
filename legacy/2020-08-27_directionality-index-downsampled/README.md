# Directionality index with downsampled contact matrices

## Motivation

Given the results in [`../2020-08-26_TAD-qc-via-downsample/`](../2020-08-26_TAD-qc-via-downsample/), despite the conclusions from multiple review papers stating that TopDom is robust to sequencing depth, we find that TAD calls between the benign and tumour samples are still affected by this feature of the data.

To compare these cases in a way that removes this confounding variable, we will downsample each file to the same number of contacts, then compare them.

## Methods

To simplify this process, we will use `cooltools diamond-insulation`, which makes use of the directionality index.

## Results

### Establishing ranges of contacts

The number of contacts in each cooler file is as follows:

| Sample Group | Min | Median | Max |
| ------------ | --- | ------ | --- |
| Primary benign | 310397759 | 385398277 | 394684060 |
| Primary tumour | 448523558 | 552585964 | 773031656 |
| Cell line benign | 207550760 | 220931955 | 234313150 |
| Cell line tumour | 124416168 | 170310595 | 302403032 |
| Cell line non-prostate | | |

The comparisons we're interested in are:

| Comparison | Shared min contacts (x10^6) |
| ---------- | --------------------------- |
| primary benign vs primary tumour | 300 |
| primary tumour vs cell line tumour | 120 |
| primary benign vs cell line benign | 200 |
| primary prostate vs cell line prostate | 120 |
| cell line prostate vs cell line non-prostate | 120 |

