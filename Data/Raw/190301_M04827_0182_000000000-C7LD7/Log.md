# 2019-03-08

I tried doing some preprocessing with Juicer, but found it to be the absolute worst piece of bioinformatics software I've ever seen.
Moreover, it does not split reads according to putative ligation sites before mapping, which will caus problems for every read that spans a ligation site.

To improve upon this, I'm switching to [HiCUP](http://www.bioinformatics.babraham.ac.uk/projects/hicup/read_the_docs/html/index.html).

## HiCUP setup

It installs via Anaconda, but it's not up to date, so I'm installing it via the download page.
See `HiC.sh` for install details.

### Creating digested reference

I've created a softlink to Bowtie2's hg19 index, and am using this for alignment.
The rule for creating the digest FASTA file is in `../../Processed/HiCUP_hg19_digested/Snakefile'.

### Processing

Instead of doing each step manually, the `hicup` command performs all steps, sequentially.
I've configured this in the `Snakefile`, and it seems to be working properly.
All the samples have completed preprocessing in a couple hours.

The duplication rates are quite high (48-77%), so this is something we'll need to be aware of beforehand.
