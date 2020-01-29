{
    FS=OFS="\t"
}
{
    # only print non-header lines
    if (NR > 1) {
        # if INV, DEL, INS, DUP
        if ($12 != "BND") {
            # chrom, start (0-indexed), end, ID
            print $1, $2 - 1, $4, $5
        # if TRA
        } else {
            # chrom_from, start_from, start_from + size, ID
            # chrom_to, start_to
            print $1, $2 - 1, $2 + $6, $5"\n"$3, $4 - 1, $4 - 1, $5
        }
    }
}