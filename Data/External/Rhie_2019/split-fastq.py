from os.path import basename
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from gzip import open as gzopen
import argparse


def split_read(read, pos=len(read) // 2):
    """
    Split read in half to produce 2 reads

    Parameters
    ----------
    read_id : Bio.SeqRecord.SeqRecord
        Read to be split
    """
    return {
        # R1/R2 as keys
        "R" + str(i + 1): SeqRecord(
            # first and second halves of sequence string
            seq=Seq(str(read.seq[(pos * i):(pos * (i + 1))]), IUPAC.ambiguous_dna),
            id=read.id,
            description="",
            letter_annotations={"phred_quality": read.letter_annotations["phred_quality"][(pos * i):(pos * (i + 1))]}
        ) for i in range(2)
    }

def split_fastq(fastq, prefix="split"):
    """
    Extract flowcell and other metadata from a FASTQ

    Parameters
    ----------
    fastq : str
        Path to input FASTQ file
    prefix : str
        Prefix for output FASTQ files (`_R1.fastq.gz` and `_R2.fastq.gz` are appended)
    """
    # check for file format
    if fastq.endswith(".gz"):
        fq_handle = gzopen(fastq, "rt")
    else:
        fq_handle = open(fastq, "r")
    out_handles = {m: gzopen(prefix + "_" + m + ".fastq.gz", "wt") for i in ["R1", "R2"]}
    # load reads for random access
    records = SeqIO.parse(fq_handle, "fastq")
    for read in chunked_records:
        # split read into each mate and fix header information
        mates = split_read(read)
        # write to each output file
        for m in mates:
            out_handles[m].write(mates[m].format("fastq"))
    fq_handle.close()
    for m in mates:
        out_handle[m].close()

if __name__ == "__main__":
    # parse command line arguments
    PARSER = argparse.ArgumentParser(description="Split a FASTQ into 2 FASTQs by splitting each read in half")
    PARSER.add_argument(
        "fq",
        type=str,
        help="Input FASTQ to split"
    )
    PARSER.add_argument(
        "-p",
        "--prefix",
        type=str,
        help="Prefix for output FASTQ files (`_R1.fastq.gz` and `_R2.fastq.gz` are appended)",
        default="split"
    )
    ARGS = PARSER.parse_args()
    # run splitting function
    split_fastq(ARGS.fq, ARGS.prefix)
