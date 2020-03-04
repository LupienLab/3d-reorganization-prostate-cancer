# this assembly method was grabbed from the higlass documentation
# https://docs.higlass.io/data_preparation.html#gene-annotation-tracks


# 1. Set assembly and species ID
# ==============================================================================
echo "Setting assembly"
ASSEMBLY="hg38"
TAXID="9606"

# Download NCBI genbank data
DATADIR="Data"
mkdir -p $DATADIR/genbank
#wget -N -P $DATADIR ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2refseq.gz
#wget -N -P $DATADIR ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz
#wget -N -P $DATADIR ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2pubmed.gz

# 2. Download reference data from UCSC and NCBI
# ==============================================================================
# Download UCSC refGene database for assembly of interest
mkdir -p $DATADIR/$ASSEMBLY
#wget -N -P $DATADIR/$ASSEMBLY/ http://hgdownload.cse.ucsc.edu/goldenPath/$ASSEMBLY/database/refGene.txt.gz

# Filter genbank data for species of interest
echo "Filtering Genbank for human data"
zcat $DATADIR/gene2refseq.gz | grep ^${TAXID} > $DATADIR/$ASSEMBLY/gene2refseq
zcat $DATADIR/gene_info.gz | grep ^${TAXID} | sort -k 2 > $DATADIR/$ASSEMBLY/gene_info
zcat $DATADIR/gene2pubmed.gz | grep ^${TAXID} > $DATADIR/$ASSEMBLY/gene2pubmed

# Sort
# Optional: filter out unplaced and unlocalized scaffolds (which have a "_" in the chrom name)
zcat $DATADIR/$ASSEMBLY/refGene.txt.gz | awk -F $'\t' '{if (!($3 ~ /_/)) print;}' | sort -k 2 > $DATADIR/$ASSEMBLY/refGene_sorted


# 3. Get full model and citation count for each gene
# ==============================================================================
echo "Getting full model"
# Count pubmed citations
# Output: {gene_id} \t {citation_count}
cat $DATADIR/$ASSEMBLY/gene2pubmed | awk '{print $2}' | sort | uniq -c | awk '{print $2 "\t" $1}' | sort > $DATADIR/$ASSEMBLY/gene2pubmed-count

# Gene2refseq dictionary
# Output: {gene_id} \t {refseq_id}
cat $DATADIR/$ASSEMBLY/gene2refseq | awk -F $'\t' '{ split($4,a,"."); if (a[1] != "-") print $2 "\t" a[1];}' | sort | uniq  > $DATADIR/$ASSEMBLY/geneid_refseqid

# Append refseq IDs to citation count table
# Output: {gene_id} \t {refseq_id} \t {citation_count}
join $DATADIR/$ASSEMBLY/geneid_refseqid $DATADIR/$ASSEMBLY/gene2pubmed-count  | sort -k2 > $DATADIR/$ASSEMBLY/geneid_refseqid_count

# Join the refseq gene model against gene IDs
# Output: {gene_id} \t {refseq_id} \t {chrom}(5) \t {strand}(6) \t {txStart}(7) \t {txEnd}(8) \t {cdsStart}(9) \t {cdsEnd}(10) \t {exonCount}(11) \t {exonStarts}(12) \t {exonEnds}(13)
join -1 2 -2 2 $DATADIR/$ASSEMBLY/geneid_refseqid_count $DATADIR/$ASSEMBLY/refGene_sorted | awk '{ print $2 "\t" $1 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t" $13 "\t" $3; }' | sort -k1 > $DATADIR/$ASSEMBLY/geneid_refGene_count

# Join citation counts against gene information
# output -> geneid \t symbol \t gene_type \t name \t citation_count
join -1 2 -2 1 -t $'\t' $DATADIR/$ASSEMBLY/gene_info $DATADIR/$ASSEMBLY/gene2pubmed-count | awk -F $'\t' '{print $1 "\t" $3 "\t" $10 "\t" $12 "\t" $16}' | sort -k1 > $DATADIR/$ASSEMBLY/gene_subinfo_citation_count

# 1: chr (chr1)
# 2: txStart (52301201) [9]
# 3: txEnd (52317145) [10]
# 4: geneName (ACVRL1)   [2]
# 5: citationCount (123) [16]
# 6: strand (+)  [8]
# 7: refseqId (NM_000020)
# 8: geneId (94) [1]
# 9: geneType (protein-coding)
# 10: geneDesc (activin A receptor type II-like 1)
# 11: cdsStart (52306258)
# 12: cdsEnd (52314677)
# 13: exonStarts (52301201,52306253,52306882,52307342,52307757,52308222,52309008,52309819,52312768,52314542,)
# 14: exonEnds (52301479,52306319,52307134,52307554,52307857,52308369,52309284,52310017,52312899,52317145,)
join -t $'\t' $DATADIR/$ASSEMBLY/gene_subinfo_citation_count $DATADIR/$ASSEMBLY/geneid_refGene_count | awk -F $'\t' '{print $7 "\t" $9 "\t" $10 "\t" $2 "\t" $16 "\t" $8 "\t" $6 "\t" $1 "\t" $3 "\t" $4 "\t" $11 "\t" $12 "\t" $14 "\t" $15}' > $DATADIR/$ASSEMBLY/geneAnnotations.bed

python exonU.py $DATADIR/$ASSEMBLY/geneAnnotations.bed > $DATADIR/$ASSEMBLY/geneAnnotationsExonUnions.bed
sort -k1,1 -V -k2,2n $DATADIR/$ASSEMBLY/geneAnnotationsExonUnions.bed > $DATADIR/$ASSEMBLY/geneAnnotationsExonUnions.sorted.bed

# 4. Create gene annotation track file that can be loaded into higlass
# ==============================================================================
# keep canonical chromosomes
awk '/^chr[0-9]{0,3}[XYM]?\t/ {print}' $DATADIR/$ASSEMBLY/geneAnnotationsExonUnions.sorted.bed > $DATADIR/$ASSEMBLY/geneAnnotationsExonUnions.sorted.canonical.bed

echo "Converting via clodius"
clodius aggregate bedfile --max-per-tile 20 --importance-column 5 --chromsizes-filename hg38.chrom.sizes --output-file $DATADIR/$ASSEMBLY/gene-annotations-${ASSEMBLY}.beddb --delimiter $'\t' $DATADIR/$ASSEMBLY/geneAnnotationsExonUnions.sorted.canonical.bed
