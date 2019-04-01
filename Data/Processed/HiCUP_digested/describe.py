'''
describe
==========

Calculate various stats about restriction sites
'''

from __future__ import division, absolute_import, print_function
import argparse
import os.path as path
import numpy as np
import pandas as pd

# ==============================================================================
# Constants
# ==============================================================================

# ==============================================================================
# Functions
# ==============================================================================

# ==============================================================================
# Main
# ==============================================================================

def main():
    '''
    Main
    '''
    df = pd.read_csv(ARGS.digest, sep='\t', skiprows=[0])
    sites = df.loc[(df['5\'_Restriction_Site'] == 'Re1') & (df['3\'_Restriction_Site'] == 'Re1'), :]
    # + 1 since start positions are 1-indexed
    print((sites.Fragment_End_Position - sites.Fragment_Start_Position + 1).describe())


if __name__ == '__main__':
    PARSER = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    PARSER.add_argument(
        'digest',
        type=str,
        help='Digest file'
    )
    ARGS = PARSER.parse_args()
    main()
