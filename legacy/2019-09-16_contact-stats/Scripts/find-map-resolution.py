'''
find-map-resolution.py
==========

Calculate the "map resolution" for a given sample.
Map resolution is defined by the matrix resolution (in bp) such that 80% of bins
have at least 1000 contacts.
See Rao _et al._, Cell, 2014 for this definition.
'''

from __future__ import division, absolute_import, print_function
import argparse
import pandas as pd
import cooler
from tqdm import tqdm

# ==============================================================================
# Main
# ==============================================================================

def main(uri):
    '''
    Main

    Parameters
    ==========
    uri : str
        (Multi-)Cooler file URI
    '''
    c = cooler.Cooler(uri)
    p = c.pixels()
    counts = pd.DataFrame()
    chrs = ["chr" + str(chrom) for chrom in list(range(1, 23)) + ["X", "Y"]]
    for chrom in tqdm(chrs):
        # fetch pixels for this chromosome
        chrom_pixels = p.fetch(chrom)
        sums = chrom_pixels.groupby("bin1_id").agg({"count": "sum"})
        counts = pd.concat([counts, sums])
    print(counts.loc[counts["count"] >= 1000, :].shape[0] / counts.shape[0])


if __name__ == '__main__':
    PARSER = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    PARSER.add_argument(
        'cool',
        type=str,
        help='Cooler file uri (possibly multi-resolution)'
    )
    ARGS = PARSER.parse_args()
    main(ARGS.cool)
