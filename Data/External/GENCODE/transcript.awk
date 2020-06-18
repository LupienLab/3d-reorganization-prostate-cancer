BEGIN{
    FS="(\\t|; )";
    OFS="\t";
}
{
    if (NR > 5 && $3 == "gene") {
        # clean the gene_id
        gsub(/gene_id "/, "", $9);
        gsub(/"/, "", $9);

        # clean the transcript_id
        gsub(/transcript_id "/, "", $10);
        gsub(/"/, "", $10);
        
        # clean the gene_name
        gsub(/gene_name "/, "", $12);
        gsub(/"/, "", $12);

        # clean the transcript_name
        gsub(/transcript_name "/, "", $14);
        gsub(/"/, "", $14);

        # chr, start, end, strand, gene_id, gene_name, transcript_id, transcript_name
        print $1, $4, $5, $7, $9, $12, $10, $14
    }
}