# Reorganization of the 3D genome pinpoints non-coding drivers of primary prostate tumors

This repository contains all the data and analysis related to [Reorganization of the 3D genome pinpoints non-coding drivers of primary prostate tumors](https://cancerres.aacrjournals.org/content/early/2021/10/11/0008-5472.CAN-21-2056).

Published version of the paper is available on [Cancer Res](https://aacrjournals.org/cancerres/article/81/23/5833/674850/Reorganization-of-the-3D-Genome-Pinpoints)
A preprint version of this article is available on [bioRxiv](https://www.biorxiv.org/content/10.1101/2021.01.05.425333v2).
A reproducible run of this work can be found on [CodeOcean](https://codeocean.com/capsule/3837124/tree).

## Usage

To download all the code, scripts, and results, use `git clone`:

```shell
git clone https://github.com/LupienLab/3d-reorganization-prostate-cancer.git
```

This does not download the raw sequencing data.
There are placeholder folders for the raw data, but the FASTQ files are available from the [European Genome-Phenome Archive](https://ega-archive.org/).

| Data Type               | EGA Accession Number |
| ----------------------- | -------------------- |
| Whole genome sequencing | EGAS00001000900      |
| RNA-seq                 | EGAS00001000900      |
| ChIP-seq (H3K27ac)      | EGAS00001002496      |
| Hi-C                    | EGAS00001005014      |

Processed data from the Hi-C sequencing data can be found on the [Gene Expression Omnibus (Accession GSE164347)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164347).

Raw Hi-C data from other studies can be found with the links and accession numbers below.

| Data                    | Repository                                                               | Accession Number |
| ----------------------- | ------------------------------------------------------------------------ | ---------------- |
| 22Rv1, RWPE1, and C4-2B | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118629)      | GSE118629        |
| H1-hESC (Rep 1)         | [4D Nucleome](https://data.4dnucleome.org/files-processed/4DNFI6HDY7WZ/) | 4DNFI6HDY7WZ     |
| H1-hESC (Rep 2)         | [4D Nucleome](https://data.4dnucleome.org/files-processed/4DNFITH978XV/) | 4DNFITH978XV     |
| HAP-1 (Rep 1)           | [4D Nucleome](https://data.4dnucleome.org/files-processed/4DNFIT64Q7A3/) | 4DNFIT64Q7A3     |
| HAP-1 (Rep 2)           | [4D Nucleome](https://data.4dnucleome.org/files-processed/4DNFINSKEZND/) | 4DNFINSKEZND     |
| GM12878 (Rep 1)         | [4D Nucleome](https://data.4dnucleome.org/files-processed/4DNFIIV4M7TF/) | 4DNFIIV4M7TF     |
| GM12878 (Rep 2)         | [4D Nucleome](https://data.4dnucleome.org/files-processed/4DNFIXVAKX9Q/) | 4DNFIXVAKX9Q     |

## Project Structure

This repository is structured as follows:

```shell
.
└── data/            # directory where all non-analysis data is stored
    ├── External/    # data from other papers, collaborators
    ├── Raw/         # raw data generated for this specific project along with pre-processing scripts and data
    └── Processed/   # data from `Raw/` that has been aggregated or processed in some other way beyond the standard raw pre-processing
└── code/
    ├── Result1/     # analysis scripts and logs for `result1`
    ├── Result2/     # analysis scripts and logs for `result2`
    └── ...
└── results/
    ├── Result1/     # results for `result1`
    ├── Result2/     # results for `result2`
    └── ...
├── README.md        # this file
└── environment.yaml # Anaconda environment YAML file for the entire project
```

To re-run any of the analyses in the `code/` folders:

1. Build and activate the `conda` environment stored in `environment.yaml`
    ```shell
    conda create --file environment.yaml -n <ENV_NAME>
    conda activate <ENV_NAME>
    ```
2. Navigate to the result directory of interest
3. Run `snakemake`

That should regenerate the entire set of results for that specific folder.
You can preview that needs to be run by running `snakemake -n` before running the analyses.
