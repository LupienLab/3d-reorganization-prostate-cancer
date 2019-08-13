from __future__ import division, print_function
import sys

import numpy as np
import h5py

import pandas as pd
import argparse

from cooler import Cooler
from cooler import util


def get_matrix_size(c, row_region, col_region):
    nrows = c.extent(row_region)[1] - c.extent(row_region)[0]
    ncols = c.extent(col_region)[1] - c.extent(col_region)[0]
    return ncols * nrows


def load_matrix(c, row_region, col_region, field, balanced, scale):
    mat = (c.matrix(balance=balanced, field=field)
            .fetch(row_region, col_region))
    if scale == 'log2':
        mat = np.log2(mat)
    elif scale == 'log10':
        mat = np.log10(mat)
    return mat

def filter_tad_range(tads, row_lo, row_hi, col_lo, col_hi):
    # ensure tads fit within row and col boundaries
    row_idx = (row_lo < tads['from.coord']) & (tads['to.coord'] < row_hi)
    col_idx = (col_lo < tads['from.coord']) & (tads['to.coord'] < col_hi)
    idx = row_idx & col_idx
    return tads.loc[idx, :]

def tads_to_patches(tads):
    from matplotlib.path import Path
    import matplotlib.patches as patches
    paths = []
    colours = []
    for i, t in tads.iterrows():
        if t.tag == 'gap':
            coords = [
                (t['from.coord'], t['from.coord']),
                (t['to.coord'], t['to.coord'])
            ]
            codes = (Path.MOVETO, Path.LINETO)
            paths.append(Path(coords, codes))
            colours.append('#bdbdbd')
        elif t.tag == 'boundary':
            coords = [
                (t['from.coord'], t['to.coord']),
                (t['to.coord'], t['from.coord'])
            ]
            codes = (Path.MOVETO, Path.LINETO)
            paths.append(Path(coords, codes))
            colours.append('#263238')
        elif t.tag == 'domain':
            coords = [
                (t['from.coord'], t['from.coord']),
                (t['from.coord'], t['to.coord']),
                (t['to.coord'], t['to.coord'])
            ]
            codes = (Path.MOVETO, Path.LINETO, Path.LINETO)
            paths.append(Path(coords, codes))
            colours.append('#1565c0')
    tad_patches = [patches.PathPatch(p, edgecolor=c, facecolor='none', lw=3) for p, c in zip(paths, colours)]
    return tad_patches

def show(cool_uri, range, balanced=True, out=None, dpi=300, scale='log10', force=False, zmin=None, zmax=None, cmap='YlOrRd', field='count', tads=None, verbose=None):
    try:
        import matplotlib as mpl
        if out is not None:
            mpl.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.path import Path
        import matplotlib.patches as patches
    except ImportError:
        print("Install matplotlib to use cooler show", file=sys.stderr)
        sys.exit(1)

    if verbose is not None:
        print('Loading contact matrix')
    c = Cooler(cool_uri)

    chromsizes = c.chromsizes
    row_region = range
    col_region = range
    row_chrom, row_lo, row_hi = util.parse_region(row_region, chromsizes)
    col_chrom, col_lo, col_hi = util.parse_region(col_region, chromsizes)

    if verbose is not None:
        print('Plotting contact matrix')
    fig, ax = plt.subplots(figsize=(11, 10))
    # ax.gcf().canvas.set_window_title('Contact matrix'.format())
    # ax.title('')
    im = ax.imshow(
        load_matrix(c, row_region, col_region, field, balanced, scale),
        interpolation='none',
        extent=[col_lo, col_hi, row_hi, row_lo],
        vmin=zmin,
        vmax=zmax,
        cmap=cmap
    )
    if tads is not None:
        if verbose is not None:
            print('Loading TADs')
        tads_data = pd.read_csv(tads, sep='\t', header=[0])
        tads_data = filter_tad_range(tads_data, row_lo, row_hi, col_lo, col_hi)
        tad_patches = tads_to_patches(tads_data)
        if verbose is not None:
            print('Plotting')
        for p in tad_patches:
            ax.add_patch(p)

    # If plotting into a file, plot and quit
    plt.ylabel('{} coordinate'.format(row_chrom))
    plt.xlabel('{} coordinate'.format(col_chrom))
    cb = fig.colorbar(im, cmap=cmap)
    cb.set_label(
        {'linear': 'relative contact frequency',
         'log2'  : 'log 2 ( relative contact frequency )',
         'log10' : 'log 10 ( relative contact frequency )'}[scale])
    fig.savefig(out, dpi=dpi)

if __name__ == '__main__':
    PARSER = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    PARSER.add_argument(
        'cool',
        type=str,
        help='Cooler file/URI to plot'
    )
    PARSER.add_argument(
        'range',
        type=str,
        help='Genomic range to plot (UCSC format)'
    )
    PARSER.add_argument(
        '-t' ,'--tads',
        type=str,
        help='TSV containing information about TADs'
    )
    PARSER.add_argument(
        '-o', '--output',
        type=str,
        help='Path to output image'
    )
    PARSER.add_argument(
        '-v',
        help='Verbose output',
        dest='verbose',
        action='count'
    )
    ARGS = PARSER.parse_args()

    # plot region
    show(ARGS.cool, ARGS.range, balanced=True, out=ARGS.output, dpi=300, tads=ARGS.tads, verbose=ARGS.verbose)


