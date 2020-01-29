# Summary

This folder contains the filtered BAMs and other associated files for LNCaP CTCF ChIP-seq from ENCODE:

* [ENCSR315NAC](https://www.encodeproject.org/experiments/ENCSR315NAC/)
* [ENCSR000DMF](https://www.encodeproject.org/experiments/ENCSR000DMF/)

These are the only two experiments (as of 2019-08-26) that have CTCF ChIP-seq data for LNCaP cells on GRCh38, without drug treatments.
Due to the ENCODE download links changing, I manually downloaded these before running `snakemake`.

| Experiment Accession                                                  | File Type           | Source                                                                                       |
| --------------------------------------------------------------------- | ------------------- | -------------------------------------------------------------------------------------------- |
| [ENCSR315NAC](https://www.encodeproject.org/experiments/ENCSR315NAC/) | Filtered BAM (Rep1) | [ENCFF500VJD](https://www.encodeproject.org/files/ENCFF500VJD/@@download/ENCFF500VJD.bam)    |
|                                                                       | Filtered BAM (Rep2) | [ENCFF117GZZ](https://www.encodeproject.org/files/ENCFF117GZZ/@@download/ENCFF117GZZ.bam)    |
|                                                                       | Conservative peaks  | [ENCFF292WTB](https://www.encodeproject.org/files/ENCFF292WTB/@@download/ENCFF292WTB.bigBed) |
|                                                                       | Optimal peaks       | [ENCFF700QXT](https://www.encodeproject.org/files/ENCFF700QXT/@@download/ENCFF700QXT.bigBed) |
| [ENCSR000DMF](https://www.encodeproject.org/experiments/ENCSR000DMF/) | Filtered BAM (Rep1) | [ENCFF925XUX](https://www.encodeproject.org/files/ENCFF925XUX/@@download/ENCFF925XUX.bam)    |
|                                                                       | Filtered BAM (Rep2) | [ENCFF265TUU](https://www.encodeproject.org/files/ENCFF265TUU/@@download/ENCFF265TUU.bam)    |
|                                                                       | Conservative peaks  | [ENCFF058ZMC](https://www.encodeproject.org/files/ENCFF058ZMC/@@download/ENCFF058ZMC.bigBed) |
|                                                                       | Optimal peaks       | [ENCFF957KCI](https://www.encodeproject.org/files/ENCFF957KCI/@@download/ENCFF957KCI.bigBed) |
