"""
split-sra
==========

Split FASTQs containing paired-end sequencing data where the mates are concatenated on the same line.
This is not the same as an interspersed FASTQ
"""

from __future__ import division, absolute_import, print_function
import argparse
import os.path as path
from Bio import SeqIO
from gzip import open as gzopen
from tqdm import tqdm

# ==============================================================================
# Constants
# ==============================================================================

# ==============================================================================
# Main
# ==============================================================================


def main(fq, prefix, l):
    """
    Main
    """
    compressed = False
    # extract input file name
    if prefix is None:
        prefix = path.basename(fq)
        if prefix.endswith(".gz"):
            compressed = True
            prefix = prefix.rstrip(".gz")
        if prefix.endswith(".fq"):
            prefix = prefix.rstrip(".fq")
        elif prefix.endswith(".fastq"):
            prefix = prefix.rstrip(".fastq")
    print("Writing to", prefix + "_R1.fastq.gz and", prefix + "_R2.fastq.gz")
    # create a file handle for the input FASTQ
    if compressed:
        fh_in = gzopen(fq, "rt")
    else:
        fh_in = open(fq, "r")
    # create output R1 and R2 files
    # fh_out = {i: gzopen(prefix + "_R" + str(i) + ".fastq.gz", "wt") for i in [1, 2]}
    records = SeqIO.parse(fh_in, "fastq")
    total = 0
    for r in tqdm(records):
        print(r)
        # total += 1
        # if total > 10:
        #     break
        # # remove "length=NNN" suffix from read ID
        # new_r = {
        #     1: SeqIO.SeqRecord(seq=r.seq[0:l], id=r.id, name=r.name),
        #     2: SeqIO.SeqRecord(seq=r.seq[l:], id=r.id, name=r.name),
        # }
        # print(new_r)
    #     for i in [1, 2]:
    #         SeqIO.write(new_r[i], fh_out[i], "fastq")
    # for i in [1, 2]:
    #     fh_out[i].close()
    fh_in.close()


if __name__ == "__main__":
    PARSER = argparse.ArgumentParser()
    PARSER.add_argument("fq", type=str, help="Input FASTQ to split")
    PARSER.add_argument(
        "length", type=int, help="Sequencing length of each mate.",
    )
    PARSER.add_argument(
        "-o",
        "--prefix",
        type=str,
        help="Prefix for output FASTQs. Default is to append `_R1` and `_R2` to the input name",
        default=None,
    )
    ARGS = PARSER.parse_args()
    main(ARGS.fq, ARGS.prefix, ARGS.length)
