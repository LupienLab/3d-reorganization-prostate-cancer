# this assembly method was grabbed from the higlass documentation
# https://docs.higlass.io/data_preparation.html#gene-annotation-tracks


# 1. Set assembly and species ID
# ==============================================================================
echo "Setting assembly"
ASSEMBLY="hg38"
TAXID="9606"

# 2. Download reference data from UCSC and NCBI
# ==============================================================================
# Download UCSC refGene database for assembly of interest
mkdir -p $DATADIR/$ASSEMBLY
wget -N -P $DATADIR/$ASSEMBLY http://hgdownload.cse.ucsc.edu/goldenPath/$ASSEMBLY/database/refGene.txt.gz

# Download NCBI genbank data
DATADIR="Data"
wget -N -P $DATADIR ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2refseq.gz

# Filter genbank data for species of interest
echo "Filtering Genbank for human data"
zcat $DATADIR/gene2refseq.gz | grep ^${TAXID} > $DATADIR/$ASSEMBLY/gene2refseq

# Sort
# Optional: filter out unplaced and unlocalized scaffolds (which have a "_" in the chrom name)
# 1: chr            (chr1)
# 2: txStart        (52301201) [9]
# 3: txEnd          (52317145) [10]
# 4: geneName       (ACVRL1)   [2]
# 5: citationCount  (123) [16]
# 6: strand         (+)  [8]
# 7: refseqId       (NM_000020)
# 8: geneId         (94)
# 9: geneType       (protein-coding)
# 10: geneDesc      (activin A receptor type II-like 1)
# 11: cdsStart      (52306258)
# 12: cdsEnd        (52314677)
# 13: exonStarts    (52301201,52306253,52306882,52307342,52307757,52308222,52309008,52309819,52312768,52314542,)
# 14: exonEnds      (52301479,52306319,52307134,52307554,52307857,52308369,52309284,52310017,52312899,52317145,)
zcat $DATADIR/$ASSEMBLY/refGene.txt.gz | awk '{FS=OFS="\t"}{if (!($3 ~ /_/)) print $3, $5, $6, $13, ".", $4, $2, $2, ".", ".", $7, $8, $10, $11;}' | sort -k1,1 -V -k2,2n > $DATADIR/$ASSEMBLY/refGene_sorted

# 3. Get full model and citation count for each gene
# ==============================================================================
echo "Linking exons together"

python exonU.py $DATADIR/$ASSEMBLY/refGene_sorted > $DATADIR/$ASSEMBLY/geneAnnotationsExonUnions.bed
sort -k1,1 -V -k2,2n $DATADIR/$ASSEMBLY/geneAnnotationsExonUnions.bed > $DATADIR/$ASSEMBLY/geneAnnotationsExonUnions.sorted.bed

# 4. Create gene annotation track file that can be loaded into higlass
# ==============================================================================

echo "Converting via clodius"
clodius aggregate bedfile \
    --max-per-tile 10 \
    --chromsizes-filename hg38.chrom.sizes \
    --output-file $DATADIR/$ASSEMBLY/gene-annotations-${ASSEMBLY}.beddb \
    --delimiter $'\t' \
    $DATADIR/$ASSEMBLY/geneAnnotationsExonUnions.sorted.bed
