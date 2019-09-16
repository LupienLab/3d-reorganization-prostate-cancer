# Summary

This folder contains publicly available clinical data for the CPC-GENE cohort.

| File                                                  | Source                                                                                  | Description                                                                 |
| ----------------------------------------------------- | --------------------------------------------------------------------------------------- | --------------------------------------------------------------------------- |
| `CPCG_Proteomics_MassSpec_Kislinger_PMID30889379.tsv` | [Sinha, Huang _et al._, Cancer Cell, 2019](https://doi.org/10.1016/j.ccell.2019.02.005) | Table S2. Proteins detected by mass spec in CPC-GENE patients.              |
| `CPC-GENE_Fraser-2017-Nature_Clinical-Data.tsv`       | [Fraser _et al._, Nature, 2017](https://doi.org/10.1038/nature20788)                    | Table S1. Clinical information from patients in the CPC-GENE cohort.        |
| `CPC-GENE_RNAseq_rsem_gene_FPKM.tsv`                  | [EGAD00001004424](https://www.ebi.ac.uk/ega/datasets/EGAD00001004424)                   | RNA-seq for 144 CPC-GENE samples in FPKM format.                            |
| `CPC-GENE_RNAseq_rsem_gene_FPKM_13-Low-C-only.tsv`    | [EGAD00001004424](https://www.ebi.ac.uk/ega/datasets/EGAD00001004424)                   | RNA-seq for the 13 CPC-GENE samples for which we have Low-C in FPKM format. |

For the protein data, the columns are as follows:

| Column Label  | Definition                                           |
| ------------- | ---------------------------------------------------- |
| Protein Group | Groups of protein defined by UniProt accession id    |
| Gene names    | Gene symbol                                          |
| Protein name  | Full protein name, as per UniProt                    |
| CPCG...       | Quantification of protein groups in 76 tumor samples |
