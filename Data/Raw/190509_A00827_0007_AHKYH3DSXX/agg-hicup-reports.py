'''
agg-hicup-reports.py
==========

HiCUP produces a summary report for the entire pipeline only when `hicup` is
run, not when you run the individual commands. This scripts uses the aggregated
tables I've created in `{hicup_tool}_summary.tsv` to create these summary
reports for integration with MultiQC
'''

from __future__ import division, absolute_import, print_function
import argparse
import os.path as path
import numpy as np
import pandas as pd

# ==============================================================================
# Constants
# ==============================================================================
AGG_COLS = [
    'File', 'Total_Reads_1', 'Total_Reads_2',
    'Not_Truncated_Reads_1', 'Not_Truncated_Reads_2',
    'Truncated_Read_1', 'Truncated_Read_2',
    'Average_Length_Truncated_1', 'Average_Length_Truncated_2',
    'Too_Short_To_Map_Read_1', 'Too_Short_To_Map_Read_2',
    'Unique_Alignments_Read_1', 'Unique_Alignments_Read_2',
    'Multiple_Alignments_Read_1', 'Multiple_Alignments_Read_2',
    'Failed_To_Align_Read_1', 'Failed_To_Align_Read_2',
    'Paired_Read_1', 'Paired_Read_2',
    'Valid_Pairs', 'Invalid_Pairs', 'Same_Circularised', 'Same_Dangling_Ends',
    'Same_Fragment_Internal', 'Re_Ligation', 'Contiguous_Sequence',
    'Wrong_Size',
    'Deduplication_Read_Pairs_Uniques', 'Deduplication_Cis_Close_Uniques',
    'Deduplication_Cis_Far_Uniques', 'Deduplication_Trans_Uniques',
    'Percentage_Mapped', 'Percentage_Valid', 'Percentage_Uniques',
    'Percentage_Unique_Trans', 'Percentage_Ditags_Passed_Through_HiCUP'
]

# ==============================================================================
# Functions
# ==============================================================================

# ==============================================================================
# Main
# ==============================================================================


def main(cfg):
    '''
    Main
    '''
    # load sample metadata
    metadata = pd.read_csv(cfg, sep='\t', index_col=False)
    samples = metadata.Sample.tolist()
    # load summary files
    truncation = pd.read_csv(
        path.join('Reports', 'hicup_truncater_summary.tsv'),
        sep='\t', index_col=False
    )
    mapping = pd.read_csv(
        path.join('Reports', 'hicup_mapper_summary.tsv'),
        sep='\t', index_col=False
    )
    filtering = pd.read_csv(
        path.join('Reports', 'hicup_filter_summary.tsv'),
        sep='\t', index_col=False
    )
    dedup = pd.read_csv(
        path.join('Reports', 'hicup_deduplicator_summary.tsv'),
        sep='\t', index_col=False
    )
    # aggregated table
    agg_data = pd.DataFrame(columns=AGG_COLS)
    agg_data.loc[:, 'File'] = samples
    # iterate through samples and extract meaningful results
    for i, r in agg_data.iterrows():
        s = r.at['File']
        # extract sample-specific data from aggregate tables
        sample_trunc = truncation.loc[truncation.File.str.startswith(s)]
        sample_map = mapping.loc[mapping.File.str.startswith(s)]
        sample_filter = filtering.loc[filtering.File.str.startswith(s)]
        sample_dedup = dedup.loc[dedup.File.str.startswith(s)]
        # assign information to agg_data
        # what are the column names in the truncation file
        individual_trunc_cols = [
            'Total_Reads_Processed', 'Truncated', 'Not_truncated',
            'Average_length_truncated_sequence'
        ]
        # what are the column names in the mapping file
        individual_map_cols = [
            'Reads_too_short_to_map',
            'Unique_alignments', 'Multiple_alignments', 'Failed_to_align',
            'Paired'
        ]
        # what are the relevant names in the aggregated file
        pair_trunc_cols = [
            'Total_Reads', 'Truncated_Read', 'Not_Truncated_Reads',
            'Average_Length_Truncated'
        ]
        pair_map_cols = [
            'Too_Short_To_Map_Read', 'Unique_Alignments_Read',
            'Multiple_Alignments_Read', 'Failed_To_Align_Read', 'Paired_Read'
        ]
        # bring column data together for truncation data
        for indiv, pair in zip(individual_trunc_cols, pair_trunc_cols):
            r[pair + '_1'] = sample_trunc.loc[:, indiv].values[0]
            r[pair + '_2'] = sample_trunc.loc[:, indiv].values[1]
        # bring column data together for mapping data
        for indiv, pair in zip(individual_map_cols, pair_map_cols):
            r[pair + '_1'] = sample_map.loc[:, indiv].values[0]
            r[pair + '_2'] = sample_map.loc[:, indiv].values[1]
        # aggregate filtering data
        individual_filt_cols = [
            'Valid_pairs', 'Invalid_pairs', 'Same_circularised', 'Same_dangling_ends', 'Same_internal', 'Re-ligation', 'Contiguous_sequence', 'Wrong_size'
        ]
        agg_filt_cols = [
            'Valid_Pairs', 'Invalid_Pairs', 'Same_Circularised', 'Same_Dangling_Ends', 'Same_Internal', 'Re_Ligation', 'Contiguous_Sequence', 'Wrong_Size'
        ]
    return agg_data


if __name__ == '__main__':
    PARSER = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    PARSER.add_argument(
        'cfg',
        type=str,
        help='Configuration metadata file',
        default='config.tsv'
    )
    ARGS = PARSER.parse_args()
    main(ARGS.cfg)
