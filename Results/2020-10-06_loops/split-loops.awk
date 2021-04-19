BEGIN{
    FS=OFS="\t";
    d="Loops/by-sample";
}
{
    if (NR > 1) {
        # column order:
        # 1 = Sample_ID
        # 2 = Type (Malignant or Benign)
        # 3 = chr_from
        # 4 = start_from
        # 5 = end_from
        # 6 = chr_to
        # 7 = start_to
        # 8 = end_to
        # 9 = loop_ID
        # 10 = anchor_ID_x
        # 11 = anchor_ID_y
        # 12 = detection_scale
        # 13 = fdr
        print $3, $4, $5, $6, $7, $8, $11, $9, $10, $12, $13 > d"/"$1".loops.bedpe"
    }
}
