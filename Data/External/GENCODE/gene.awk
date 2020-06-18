BEGIN{
    FS="(\\t|; )";
    OFS="\t";
}
{
    if (NR > 5 && $3 == "gene") {
        gsub(/gene_id "/, "", $9);
        gsub(/"/, "", $9);
        gsub(/gene_name "/, "", $11);
        gsub(/"/, "", $11);
        print $1, $4, $5, $7, $9, $11
    }
}