"""
bf2vcf.py
==========

Convert hic_breakfinder output files to the Variant Call Format (VCF)
"""

from __future__ import division, absolute_import, print_function
import argparse
import pandas as pd
from datetime import date
import numpy as np

# ==============================================================================
# Constants
# ==============================================================================
VCF_HEADER = """##fileformat=VCFv4.3
##fileDate={date}
##contig=<ID=chr1,length=248956422>
##contig=<ID=chr2,length=242193529>
##contig=<ID=chr3,length=198295559>
##contig=<ID=chr4,length=190214555>
##contig=<ID=chr5,length=181538259>
##contig=<ID=chr6,length=170805979>
##contig=<ID=chr7,length=159345973>
##contig=<ID=chr8,length=145138636>
##contig=<ID=chr9,length=138394717>
##contig=<ID=chr10,length=133797422>
##contig=<ID=chr11,length=135086622>
##contig=<ID=chr12,length=133275309>
##contig=<ID=chr13,length=114364328>
##contig=<ID=chr14,length=107043718>
##contig=<ID=chr15,length=101991189>
##contig=<ID=chr16,length=90338345>
##contig=<ID=chr17,length=83257441>
##contig=<ID=chr18,length=80373285>
##contig=<ID=chr19,length=58617616>
##contig=<ID=chr20,length=64444167>
##contig=<ID=chr21,length=46709983>
##contig=<ID=chr22,length=50818468>
##contig=<ID=chrX,length=156040895>
##contig=<ID=chrY,length=57227415>
##contig=<ID=chrM,length=16569>
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=TRA,Description="Translocation">
##ALT=<ID=INS,Description="Insertion">
##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for END coordinate in case of a translocation">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of the SV">
##INFO=<ID=SVRES,Number=1,Type=Integer,Description="Resolution of SV call in kbp">
##INFO=<ID=LO,Number=1,Type=Float,Description="Log10 odds ratio for SV call">
##INFO=<ID=CT,Number=1,Type=String,Description="Paired-end signature induced connection type">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
"""

# ==============================================================================
# Main
# ==============================================================================


def main(infile, outfile):
    """
    Main
    """
    bf_calls = pd.read_csv(
        infile,
        sep="\t",
        index_col=False,
        names=[
            "logodds",
            "chr_from",
            "start_from",
            "end_from",
            "strand_from",
            "chr_to",
            "start_to",
            "end_to",
            "strand_to",
            "resolution",
        ],
    )
    fout = open(outfile, "w")
    fout.write(VCF_HEADER.format(date=date.today().isoformat()))
    # convert from natural log (as called in breakfinder) to log10(odds ratio)
    qual_conv = np.log10(np.exp(1))
    for row in bf_calls.itertuples():
        # annotate as translocation if on different chromosomes
        if row.chr_from != row.chr_to:
            alt = "<TRA>"
        else:
            alt = "."
        # remove "kb" from Breakfinder resolution column
        if row.resolution.endswith("kb"):
            res = row.resolution.rstrip("kb")
        elif row.resolution.endswith("Mb"):
            res = 1000 * int(row.resolution.rstrip("Mb"))
        output = [
            # chrom
            row.chr_from,
            # pos
            row.start_from,
            # id
            ".",
            # ref
            "N",
            # alt
            alt,
            # qual
            ".",
            # filter
            "PASS",
            # info
            "CHR2={chrom};END={end};SVRES={res};SVLEN={length};LO={lo}".format(
                chrom=row.chr_to,
                end=row.start_to,
                res=res,
                length=row.end_from - row.start_from,
                lo=np.round(qual_conv * row.logodds, 3),
            ),
        ]
        fout.write("\t".join([str(o) for o in output]) + "\n")
    fout.close()


if __name__ == "__main__":
    PARSER = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    PARSER.add_argument(
        "bf", type=str, help="`hic_breakfinder` output file to be converted"
    )
    PARSER.add_argument(
        "-o", "--output", type=str, help="Output file name", default="breaks.vcf"
    )
    ARGS = PARSER.parse_args()
    main(ARGS.bf, ARGS.output)
