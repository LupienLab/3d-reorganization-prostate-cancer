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
    'File',
    'Total_Reads_1',
    'Total_Reads_2',
    'Not_Truncated_Reads_1',
    'Not_Truncated_Reads_2',
    'Truncated_Read_1',
    'Truncated_Read_2',
    'Average_Length_Truncated_1',
    'Average_Length_Truncated_2',
    'Too_Short_To_Map_Read_1',
    'Too_Short_To_Map_Read_2',
    'Unique_Alignments_Read_1',
    'Unique_Alignments_Read_2',
    'Multiple_Alignments_Read_1',
    'Multiple_Alignments_Read_2',
    'Failed_To_Align_Read_1',
    'Failed_To_Align_Read_2',
    'Paired_Read_1',
    'Paired_Read_2',
    'Valid_Pairs',
    'Invalid_Pairs',
    'Same_Circularised',
    'Same_Dangling_Ends',
    'Same_Fragment_Internal',
    'Re_Ligation',
    'Contiguous_Sequence',
    'Wrong_Size',
    'Deduplication_Read_Pairs_Uniques',
    'Deduplication_Cis_Close_Uniques',
    'Deduplication_Cis_Far_Uniques',
    'Deduplication_Trans_Uniques',
    'Percentage_Mapped',
    'Percentage_Valid',
    'Percentage_Uniques',
    'Percentage_Unique_Trans',
    'Percentage_Ditags_Passed_Through_HiCUP'
]

BATCHES = ['190509_A00827_0007_AHKYH3DSXX', '190605_A00827_0009_BHLGJ3DSXX']

# ==============================================================================
# Functions
# ==============================================================================


def read_batch_tsv(batch):
    '''
    Simple wrapper for pd.read_csv

    Parameters
    ----------
    batch : str
        Sequencing batch data folder
    '''
    raw_data_dir = path.join('..', '..', 'Raw')
    report_dir = 'Reports'
    summary_name = 'HiCUP_summary_report.tsv'
    return pd.read_csv(
        path.join(raw_data_dir, batch, report_dir, summary_name),
        sep='\t', index_col=False
    )

# ==============================================================================
# Main
# ==============================================================================


def main(cfg, outdir):
    '''
    Main
    '''
    # load sample metadata
    metadata = pd.read_csv(cfg, sep='\t', index_col=False)
    samples = metadata.Sample.tolist()
    frames = [read_batch_tsv(b) for b in BATCHES]
    agg_data = pd.concat(frames)
    # add sample column to table for aggregation
    agg_data['Sample'] = None
    for s in samples:
        agg_data.loc[agg_data.File.str.startswith(s), 'Sample'] = s
    # aggregate by sample, add all count columns
    full_data = agg_data.loc[:, ~(agg_data.columns.str.startswith(
        'Percentage') | agg_data.columns.str.startswith('Average'))].groupby('Sample').agg('sum')

    # calculate percentages (ignoring the average lengths, not necessary)
    full_data['Percentage_Mapped'] = np.round(
        100 * (full_data['Paired_Read_1'] + full_data['Paired_Read_2']) / (
            full_data['Total_Reads_1'] + full_data['Total_Reads_2']),
        decimals=2
    )
    full_data['Percentage_Valid'] = np.round(
        100 * full_data['Valid_Pairs'] /
        (full_data['Valid_Pairs'] + full_data['Invalid_Pairs']),
        decimals=2
    )
    full_data['Percentage_Uniques'] = np.round(
        100 * full_data['Deduplication_Read_Pairs_Uniques'] /
        full_data['Valid_Pairs'],
        decimals=2
    )
    full_data['Percentage_Unique_Trans'] = np.round(
        100 * full_data['Deduplication_Trans_Uniques'] /
        full_data['Deduplication_Read_Pairs_Uniques'],
        decimals=2
    )
    full_data['Percentage_Ditags_Passed_Through_HiCUP'] = np.round(
        200 * full_data['Deduplication_Read_Pairs_Uniques'] /
        (full_data['Total_Reads_1'] + full_data['Total_Reads_2']),
        decimals=2
    )
    # rename `Sample` to `File` to match HiCUP report format
    full_data.reset_index(drop=False, inplace=True)
    full_data.columns.values[0] = 'File'
    # save data to output directory
    full_data.to_csv(
        path.join(outdir, 'HiCUP_summary_report.tsv'),
        index=False, sep='\t', header=True
    )
    for i, r in full_data.iterrows():
        full_data.loc[[i]].to_csv(
            path.join(outdir, 'HiCUP_summary_report_{}.txt'.format(r.File)),
            index=False, sep='\t', header=True
        )


if __name__ == '__main__':
    PARSER = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    PARSER.add_argument(
        'cfg',
        type=str,
        help='Configuration metadata file'
    )
    PARSER.add_argument(
        '-o', '--outdir',
        type=str,
        help='Output directory for summary files',
        default='.'
    )
    ARGS = PARSER.parse_args()
    main(ARGS.cfg, ARGS.outdir)
