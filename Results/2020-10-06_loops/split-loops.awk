BEGIN{
    FS=OFS="\t";
    d="Loops/by-sample";
}
{
    if (NR > 1) {
        print $4, $5, $6, $7, $8, $9, $1, $10, $11, $12, $13 > d"/"$2".loops.bedpe"
    }
}
