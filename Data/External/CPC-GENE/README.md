# Summary

This folder contains publicly available clinical data for the CPC-GENE cohort.

## Datasets

| File | Source | Description |
| ---- | ------ | ----------- |
| `CPCG_Proteomics_MassSpec_Kislinger_PMID30889379.tsv` | [Sinha, Huang _et al._, Cancer Cell, 2019](https://doi.org/10.1016/j.ccell.2019.02.005) | Table S2. Proteins detected by mass spec in CPC-GENE patients. |
| `CPC-GENE_Fraser-2017-Nature_Clinical-Data.tsv`       | [Fraser _et al._, Nature, 2017](https://doi.org/10.1038/nature20788)                    | Table S1. Clinical information from patients in the CPC-GENE cohort. |
| `CPC-GENE_RNAseq_rsem_gene_FPKM.tsv`                  | [Chen, Huang, Xu, Livingstone, Soares, _et al._, Cell, 2019](https:/doi.org/10.1016/j.cell.2019.01.025) | RNA-seq for 144 CPC-GENE samples in FPKM units from RSEM. Raw data from [EGAD00001004424](https://www.ebi.ac.uk/ega/datasets/EGAD00001004424). |
| `CPC-GENE_RNAseq_rsem_gene_FPKM_13-Low-C-only.tsv`                  | [Chen, Huang, Xu, Livingstone, Soares, _et al._, Cell, 2019](https:/doi.org/10.1016/j.cell.2019.01.025) | Same as above, but only for the 13 samples for which we have performed the Low-C protocol on. |
| `CPC-GENE_Fraser-2017_SVs-Mbp-bins.tsv`               | [Fraser _et al._, Nature, 2017](https://doi.org/10.1038/nature20788)                    | SVs detected in each patient, grouped into Mbp bins. Inversions = 1, deletions = 2, inter-chromosomal translocations = 3. `INV`, `DEL`, and `CTX` columns list percentage of patients with this type of structural variant in this bin. `Counts` column lists the total number of patients where some SV was detected in that bin. Coordinates in hg19. |
| `CPC-GENE_Fraser-2017_SVs-Mbp-bins_13-Low-C-only.tsv` | [Fraser _et al._, Nature, 2017](https://doi.org/10.1038/nature20788)                    | Same as above, but only for the 13 samples for which we have performed the Low-C protocol on. |

For the protein data, the columns are as follows:

| Column Label  | Definition                                           |
| ------------- | ---------------------------------------------------- |
| Protein Group | Groups of protein defined by UniProt accession id    |
| Gene names    | Gene symbol                                          |
| Protein name  | Full protein name, as per UniProt                    |
| CPCG...       | Quantification of protein groups in 76 tumor samples |

All genome coordinates from these data are listed in hg19, unless otherwise specified.
The rest of this work is in hg38, so we've made sure to convert or liftover to hg38 coordinates where needed.
