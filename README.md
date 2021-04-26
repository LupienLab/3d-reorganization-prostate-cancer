# Cis-Regulatory Element Hijacking Overshadows Higher-Order Topological Changes in Prostate Cancer

This repository contains all data and analysis related to 

Hi-C data analysis in primary prostate cancer.

## Usage

This repository is structured as follows:

```shell
.
└── data/            # directory where all non-analysis data is stored
    ├── External/    # data from other papers, collaborators
    ├── Raw/         # raw data generated for this specific project along with pre-processing scripts and data
    └── Processed/   # data from `Raw/` that has been aggregated or processed in some other way beyond the standard raw pre-processing
└── code/
    ├── Result1/     # analysis scripts, logs, and results for `result1`
    ├── Result2/     # analysis scripts, logs, and results for `result2`
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
