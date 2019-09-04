from __future__ import division, print_function
import numpy as np
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


def filter_ranges(ranges, row_lo, row_hi, col_lo, col_hi):
    '''
    Filter genomic ranges based on whether they are within the extent being plotted

    Parameters
    ----------
    row_lo : int
        Lower bound on rows
    row_hi : int
        Upper bound on rows
    col_lo : int
        Lower bound on cols
    col_hi : int
        Upper bound on cols
    '''
    # ensure tads fit within row and col boundaries
    row_idx = (row_lo < ranges['from.coord']) & (ranges['to.coord'] < row_hi)
    col_idx = (col_lo < ranges['from.coord']) & (ranges['to.coord'] < col_hi)
    idx = row_idx & col_idx
    return ranges.loc[idx, :]


def tads_to_patches(tads):
    from matplotlib.path import Path
    import matplotlib.patches as patches
    paths = []
    colours = []
    for i, t in tads.iterrows():
        a = t['from.coord']
        b = t['to.coord']
        if t.tag == 'gap':
            coords = [
                (a, b),
                (b, b)
            ]
            codes = (Path.MOVETO, Path.LINETO)
            paths.append(Path(coords, codes))
            colours.append('#bdbdbd')
        elif t.tag == 'boundary':
            coords = [
                (a, a),
                ((a + b) / 2, (a + b) / 2),
                (a, b),
                ((a + b) / 2, (a + b) / 2),
                (b, b)
            ]
            codes = (Path.MOVETO, Path.LINETO,
                     Path.LINETO, Path.LINETO, Path.LINETO)
            paths.append(Path(coords, codes))
            colours.append('#263238')
        elif t.tag == 'domain':
            coords = [
                (a, a),
                (b, a),
                (b, b)
            ]
            codes = (Path.MOVETO, Path.LINETO, Path.LINETO)
            paths.append(Path(coords, codes))
            colours.append('#1565c0')
    tad_patches = [patches.PathPatch(p, edgecolor=c, facecolor='none', lw=2)
                   for p, c in zip(paths, colours)]
    return tad_patches


def parse_highlights(hl_row, hl_col, row_lo, row_hi, col_lo, col_hi, chromsizes):
    '''
    Parse highlights supplied as comma-separated UCSC-formatted regions

    Parameters
    ----------
    hl_row : str
        List of rows to highlight
    hl_col : str
        List of cols to highlight
    row_lo : int
        Lower bound on rows
    row_hi : int
        Upper bound on rows
    col_lo : int
        Lower bound on cols
    col_hi : int
        Upper bound on cols
    chromsizes : dict(str, int)
        Chromosome sizes to filter by
    '''
    import matplotlib.patches as patches
    hl_patches = []
    if hl_row is not None:
        rows = hl_row.split(',')
        rows = [util.parse_region(r, chromsizes) for r in hl_row.split(',')]
        n_rows = len(rows)
    else:
        rows = None
        n_rows = 0
    if hl_col is not None:
        cols = hl_col.split(',')
        cols = [util.parse_region(c, chromsizes) for c in hl_col.split(',')]
        n_cols = len(cols)
    else:
        cols = None
        n_cols = 0
    # parse rows and columns as pairs for higlighting
    if n_rows == n_cols:
        hl_patches = [
            patches.Rectangle(
                xy=(c[1], r[1]), height=(r[2] - r[1]), width=(c[2] - c[1]),
                edgecolor='#03DAC6', facecolor='none', lw=2) for r, c in zip(rows, cols)
        ]
    else:
        # get higlights in pairs
        m = min(n_rows, n_cols)
        if m > 0:
            hl_patches = [
                patches.Rectangle(
                    xy=(c[1], r[1]), height=(r[2] - r[1]), width=(c[2] - c[1]),
                    edgecolor='#03DAC6', facecolor='none', lw=2) for r, c in zip(rows[:m], cols[:m])
            ]
        # get highlights not in pairs (only one of these next for loops will run)
        #   cols first
        if n_cols > 0:
            for c in cols[m:]:
                hl_patches.append(
                    patches.Rectangle(
                        xy=(c[1], row_lo), height=(row_hi - row_lo), width=(c[2] - c[1]),
                        edgecolor='#03DAC6', facecolor='none', lw=2)
                )
        #   then rows
        if n_rows > 0:
            for r in rows[m:]:
                hl_patches.append(
                    patches.Rectangle(
                        xy=(col_lo, r[1]), height=(r[2] - r[1]), width=(col_hi - col_lo),
                        edgecolor='#03DAC6', facecolor='none', lw=2)
                )
    return hl_patches


def show(
    cool_uri, range, range2, balanced=True, out=None, dpi=300, scale='log10',
    hl_row=None, hl_col=None,
    force=False, zmin=None, zmax=None, cmap='YlOrRd', field='count',
    tads=None, rotate=False,
    verbose=None,
    strip_text=False
):
    '''
    Plot contact matrix

    Parameters
    ----------
    cool_uri : str
        URI to Cooler file to plot
    range : str
        UCSC-compatible range to plot (chr:start-end)
    range2 : str
        UCSC-compatible range to plot (chr:start-end)
    balanced : bool
        Plot balanced or raw contact matrix
    out : str
        Output image file name
    dpi : int
        Output image resolution
    scale : str
        Plotting transform function. One of ['linear', 'log2', 'log10'].
    hl_row : str
        Comma-separated list of row regions to highlight
    hl_col : str
        Comma-separated list of col regions to highlight
    force : bool
        Force plotting despite memory conerns
    zmin : float
        Min value that colour scale covers
    zmax : float
        Max value that colour scale covers
    cmap : str
        Colour map
    field : str
        Column in Cooler file table to plot
    tads : str
        Path to TSV file containing TAD calls
    rotate : bool
        Rotate plot by 45 degrees
    verbose : int
        Verbosity level
    strip_text : bool
        Strip plot of text and axes labels
    '''
    # import matplotlib package here for faster loading
    try:
        import matplotlib as mpl
        if out is not None:
            mpl.use('Agg')
        import matplotlib.pyplot as plt
        import matplotlib.transforms as mtransforms
    except ImportError:
        raise ImportError("Install matplotlib to use `cooler show`")
    # types of scaling before plotting
    scale_types = {
        'linear': 'relative contact frequency',
        'log2': 'log 2 ( relative contact frequency )',
        'log10': 'log 10 ( relative contact frequency )'
    }
    if scale not in scale_types:
        raise ValueError("Invalid scale. Expected one of: %s" %
                         scale_types.dict_keys())

    # load contact matri
    if verbose is not None:
        print('Loading contact matrix')
    c = Cooler(cool_uri)

    # parse regions specified to plot
    chromsizes = c.chromsizes
    col_region = range
    row_region = col_region if range2 is None else range2
    row_chrom, row_lo, row_hi = util.parse_region(row_region, chromsizes)
    col_chrom, col_lo, col_hi = util.parse_region(col_region, chromsizes)

    if verbose is not None:
        print('Plotting contact matrix')
    fig, ax = plt.subplots(figsize=(5, 4))
    im = ax.imshow(
        load_matrix(c, row_region, col_region, field, balanced, scale),
        interpolation='none',
        extent=[col_lo, col_hi, row_hi, row_lo],
        vmin=zmin,
        vmax=zmax,
        cmap=cmap
    )
    # transform data images if `--rotate` option
    if rotate:
        # shift_coords = ax.transData.transform([col_hi, row_lo])
        tr = (
            mtransforms.Affine2D()
            .rotate_deg(-45)
            .scale(1 / np.sqrt(2))
            .translate(0, row_hi)
        ) + ax.transData
        # hide axes
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        # hide y axis labels and ticks
        ax.axes.get_yaxis().set_visible(False)
    else:
        tr = ax.transData
    # apply transformation to all plotted objects
    #   the contact matrix
    im.set_transform(tr)
    # add annotations for TADs if present
    if tads is not None:
        if verbose is not None:
            print('Loading TADs')
        tads_data = pd.read_csv(tads, sep='\t', header=[0])
        tads_data = filter_ranges(tads_data, row_lo, row_hi, col_lo, col_hi)
        tad_patches = tads_to_patches(tads_data)
        # add TAD calls
        for p in tad_patches:
            p.set_transform(tr)
            ax.add_patch(p)
    # add annotations for specific regions pairs if present
    if hl_row is not None or hl_col is not None:
        highlights = parse_highlights(
            hl_row, hl_col,
            row_lo, row_hi, col_lo, col_hi,
            chromsizes
        )
        # add TAD calls
        for p in highlights:
            p.set_transform(tr)
            ax.add_patch(p)

    # remove text and labels if desired
    if strip_text:
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
    # If plotting into a file, plot and quit
    plt.ylabel('{} coordinate'.format(row_chrom))
    plt.xlabel('{} coordinate'.format(col_chrom))
    cb = fig.colorbar(im, cmap=cmap)
    cb.set_label(scale_types[scale])
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
        '-r2', '--range2',
        type=str,
        help='Genomic range to plot (UCSC format)'
    )
    PARSER.add_argument(
        '--hr',
        type=str,
        help='Comma-separated list of rows to highlight (UCSC format). Pair with`--hc` to highlight specific loci.'
    )
    PARSER.add_argument(
        '--hc',
        type=str,
        help='Comma-separated list of columns to highlight (UCSC format). Pair with`--hc` to highlight specific loci.'
    )
    PARSER.add_argument(
        '-t', '--tads',
        type=str,
        help='TSV containing information about TADs'
    )
    PARSER.add_argument(
        '-r', '--rotate',
        action='store_true',
        help='Rotate contact matrix by 45 degrees in plot'
    )
    PARSER.add_argument(
        '--zmin',
        type=float,
        help='Minimal value on colour scale',
        default=None
    )
    PARSER.add_argument(
        '--zmax',
        type=float,
        help='Maximal value on colour scale',
        default=None
    )
    PARSER.add_argument(
        '--no-text',
        action='store_true',
        help='Strip plots of text and axes'
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
    show(
        cool_uri=ARGS.cool, range=ARGS.range, range2=ARGS.range2, balanced=True, out=ARGS.output,
        dpi=600,
        hl_row=ARGS.hr, hl_col=ARGS.hc,
        zmin=ARGS.zmin, zmax=ARGS.zmax,
        tads=ARGS.tads, rotate=ARGS.rotate, verbose=ARGS.verbose,
        strip_text=ARGS.no_text
    )
