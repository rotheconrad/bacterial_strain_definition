 #!/usr/bin/env python

'''fastANI shared genome fraction vs ANI correlation.

This script reads fastANI output of one vs many or all vs all and builds
a scatter plot with the mean, median and correlation. The y-axis is the
shared genome fraction and the x-axis is the ANI. The fastANI output
includes the count of bidirectional fragment mappings and the total query
fragments along with the ANI estimate. The shared genome fraction is the
bidirectional fragment mapping count / the total fragments.

mean, meadian and correlation of ANI values (x-axis) and Shared / Total
Fragments (y-axis). Returns a scatter plot of the data as a pdf file.

This script takes the following input parameters:

    * ani - all vs all ani file from fastANI (str)
    * org - organism name for the plot title (str)
    * op - An output file prefix (str)

This script returns a publication ready figure in pdf format.

This script requires python 3.6+ and the following modules:

    * matplotlib
    * numpy
    * pandas
    * seaborn
    * scipy
    * pygam - https://pygam.readthedocs.io/
        - pip install pygam
        - conda install -c conda-forge pygam
        - conda install -c conda-forge scikit-sparse nose
        - goes way faster with more RAM. More RAM, more fast.
    * datashader for density scatter plot if > 500 genome pairs.
        - pip install datashader
        - conda install datashader

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: April 2022
License :: GNU GPLv3
Copyright 2022 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse, sys, random
from collections import defaultdict
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np, pandas as pd; np.random.seed(0)
import seaborn as sns; sns.set(style="white", color_codes=True)
from scipy.stats import pearsonr as corr
from pygam import LinearGAM


"""
# Set the default sans-serif font to Helvetica
matplotlib.rcParams['font.sans-serif'] = "Helvetica"
# Set default font to sans-serif
matplotlib.rcParams['font.family'] = "sans-serif"
"""

class TickRedrawer(matplotlib.artist.Artist):
    #https://stackoverflow.com/questions/19677963/
    #matplotlib-keep-grid-lines-behind-the-graph-but-the-y-and-x-axis-above
    """Artist to redraw ticks."""

    __name__ = "ticks"

    zorder = 10

    @matplotlib.artist.allow_rasterization
    def draw(self, renderer: matplotlib.backend_bases.RendererBase) -> None:
        """Draw the ticks."""
        if not self.get_visible():
            self.stale = False
            return

        renderer.open_group(self.__name__, gid=self.get_gid())

        for axis in (self.axes.xaxis, self.axes.yaxis):
            loc_min, loc_max = axis.get_view_interval()

            for tick in axis.get_major_ticks() + axis.get_minor_ticks():
                if tick.get_visible() and loc_min <= tick.get_loc() <= loc_max:
                    for artist in (tick.tick1line, tick.tick2line):
                        artist.draw(renderer)

        renderer.close_group(self.__name__)
        self.stale = False


def subsample_gpairs(gpairs, npairs, r, xmin, xmax):
    """ subsample r gpairs between [xmin, xmax] ANI """
    X = []

    while len(X) < r:
        z = random.sample(npairs, 1)[0]
        ani = gpairs[z][1]
        if ani >= xmin and ani <= xmax: X.append(z)

    return X

def gather_subsampled_data(ani_file, xmin, xmax, r, e):
    """Reads in the tsv ALLvsAll ANI """

    print("\nReading data.")
    name_dict = defaultdict(dict)
    data_dict = {'species': [], 'xs': [], 'ys': []}

    with open(ani_file, 'r') as f:
        for l in f:
            # split each line by tabs
            X = l.rstrip().split('\t')
            # get the genome pair names
            qry_genome = X[0]
            ref_genome = X[1]
            species = qry_genome.split('/')[1]
            # total frac column. Keep larger genome as reference
            tfrac = int(X[4])
            # remove self matchs
            if qry_genome == ref_genome: continue
            # compute metrics
            ani = float(X[2])
            # get the shared genome fraction
            shared = float(X[3])
            total = float(X[4])
            ratio = shared / total
            # sort genome file names and combine
            names = [qry_genome, ref_genome]
            names.sort()
            gname = '-'.join(names)
            # record values from only one genome pair.
            # Keep larger genome as reference.
            if tfrac > name_dict[species].get(gname[0], 0):
                # add genome pair to dict
                name_dict[species][gname] = [tfrac, ani, ratio]

    for i in range(e):
        # subsample r genome pairs per species and write data to arrays
        for species, gpairs in name_dict.items():
            npairs = list(gpairs.keys())
            vals = gpairs.values()
            # i[1] = ani
            testani = [i[1] for i in vals if i[1] >= xmin and i[1] <= xmax]
            if len(testani) < 30:
                #X = npairs
                continue
            else:
                X = subsample_gpairs(gpairs, npairs, r, xmin, xmax)
                
            for g in X:
                ani = gpairs[g][1]
                ratio = gpairs[g][2]
                data_dict['species'].append(species)
                data_dict['xs'].append(ani)
                data_dict['ys'].append(ratio)


    # convert to dataframe
    df = pd.DataFrame(data_dict)
    df = df[df['xs'] <= xmax] 
    df = df[df['xs'] >= xmin]
    n = len(df)

    return df, n


def gather_data(ani_file, xmin, xmax):
    """Reads in the tsv ALLvsAll ANI """

    print("\nReading data.")
    name_dict = {}
    data_dict = {'species': [], 'xs': [], 'ys': []}

    with open(ani_file, 'r') as f:
        for l in f:
            # split each line by tabs
            X = l.rstrip().split('\t')
            # get the genome pair names
            qry_genome = X[0]
            ref_genome = X[1]
            species = qry_genome.split('/')[1]
            # total frac column. Keep larger genome as reference
            tfrac = int(X[4]) 
            # remove self matchs
            if qry_genome == ref_genome: continue
            # compute metrics
            ani = float(X[2])
            # get the shared genome fraction
            shared = float(X[3])
            total = float(X[4])
            ratio = shared / total
            # sort genome file names and combine
            names = [qry_genome, ref_genome]
            names.sort()
            gname = '-'.join(names)
            # record values from only one genome pair.
            # Keep larger genome as reference.
            if tfrac > name_dict.get(gname[0], 0):
                # add genome pair to dict
                name_dict[gname] = [tfrac, ani, ratio, species]

    # write data to arrays
    for gpairs, metrics in name_dict.items():
        ani = metrics[1]
        ratio = metrics[2]
        species = metrics[3]
        data_dict['xs'].append(ani)
        data_dict['ys'].append(ratio)
        data_dict['species'].append(species)
    # convert to dataframe
    df = pd.DataFrame(data_dict)
    df = df[df['xs'] <= xmax] 
    df = df[df['xs'] >= xmin]
    n = len(df)

    total_species = len(set(data_dict['species']))
    filtered_species = len(df['species'].unique())
    print(f'\nTotal species in file: {total_species}')
    print(f'Species between {xmin}-{xmax} ANI: {filtered_species}')

    return df, n


def gather_stats(df):
    """Computes correlation, mean, and median on df columns xs and ys """

    # Compute Pearson Correlation Coefficient
    print("\nCalculating statistics.")
    pcorr = corr(df['xs'], df['ys'])

    # Compute ANI mean and median
    ani_mean = np.mean(df['xs'])
    ani_median = np.median(df['xs'])
    frag_mean = np.mean(df['ys'])
    frag_median = np.median(df['ys'])

    # Compile dictionairy
    df_stats = {
        'pcorr': pcorr,
        'ani_mean': ani_mean,
        'ani_median': ani_median,
        'frag_mean': frag_mean,
        'frag_median': frag_median
        }

    print(f"\nANI mean: {ani_mean:.2f}\nANI median: {ani_median:.2f}")
    print(f"\nFrag mean: {frag_mean:.2f}\nFrag median: {frag_median:.2f}")

    return df_stats


def fastANI_scatter_plot(
                df, n, species, outfile, xmin, xmax, xstep, p , a, z, g, c
                ):
    """Takes the data and builds the plot"""

    # Gather Stats
    df_stats = gather_stats(df)

    stats_line = (
        f"Pearson r: {round(df_stats['pcorr'][0], 2)}\n"
        f"p value: {round(df_stats['pcorr'][1], 2)}"
        )

    # Set Colors and markers
    grid_color = '#d9d9d9'
    main_color = '#933b41'
    second_color = '#737373'
    vline_color = '#000000'
    color = '#252525'
    marker = '.' #'o'

    # build plot
    gg = sns.JointGrid(x="xs", y="ys", data=df)
    print('\nComputing KDEs for marginal plots.')
    # x margin kde plot
    sns.histplot(
            x=df["xs"],
            ax=gg.ax_marg_x,
            legend=False,
            color=color,
            stat='probability'
            )
    # y margin kde plot
    sns.histplot(
            y=df["ys"],
            ax=gg.ax_marg_y,
            legend=False,
            color=color,
            stat='probability'
            )
    # main panel scatter plot
    if z: # density scatter plot with datashader
        print('\nComputing plot densities.')
        import datashader as ds
        from datashader.mpl_ext import dsshow
        dsartist = dsshow(
                        df,
                        ds.Point("xs", "ys"),
                        ds.count(),
                        norm="log",
                        aspect="auto",
                        ax=gg.ax_joint,
                        width_scale=3.,
                        height_scale=3.
                        )
        dsartist.zorder = 2.5

    else: # regular scatter plot
        print('\nPlotting data.')
        gg.ax_joint.plot(
                df["xs"],
                df["ys"],
                marker,
                ms=p,
                alpha=a,
                color=color,
                )

    if g:
        # Trendline with pyGAM
        print('\nCalculating trendline with pyGAM.')
        X = df["xs"].to_numpy()
        X = X[:, np.newaxis]
        y = df["ys"].to_list()

        gam = LinearGAM().gridsearch(X, y)
        XX = gam.generate_X_grid(term=0, n=500)

        gg.ax_joint.plot(
                    XX,
                    gam.predict(XX),
                    color='#FCEE21',
                    linestyle='--',
                    linewidth=1.0,
                    zorder=2.8
                    )
        gg.ax_joint.plot(
                    XX,
                    gam.prediction_intervals(XX, width=0.95),
                    color='#CBCB2C',
                    linestyle='--',
                    linewidth=1.0,
                    zorder=2.8
                    )
        r2 = gam.statistics_['pseudo_r2']['explained_deviance']
        GAM_line = f"GAM Pseudo R-Squared: {r2:.4f}"
        gg.ax_joint.text(
            0.75, 0.1, GAM_line,
            fontsize=10, color=second_color,
            verticalalignment='top', horizontalalignment='right',
            transform=gg.ax_joint.transAxes
            )


    # plot title, labels, text
    species_name = ' '.join(species.split('_'))
    ptitle = f'{species_name} (n={n})'
    gg.ax_marg_x.set_title(ptitle, fontsize=18, y=1.02)

    gg.ax_joint.set_xlabel(
        'Average nucleotide identity (%)',
        fontsize=12, y=-0.02
        )
    gg.ax_joint.set_ylabel(
        'Shared genome fraction', # 'Shared / total fragments'
        fontsize=12, x=-0.02
        )
    gg.ax_joint.text(
        0.25, 0.99, stats_line,
        fontsize=10, color=second_color,
        verticalalignment='top', horizontalalignment='right',
        transform=gg.ax_joint.transAxes
        )

    # set the axis parameters / style
    hstep = xstep/10
    gg.ax_joint.set_xticks(np.arange(xmin, xmax+hstep, xstep))
    gg.ax_joint.set_xlim(left=xmin-hstep, right=xmax+hstep)

    gg.ax_joint.set_yticks(np.arange(0.6, 1.1, 0.1))
    gg.ax_joint.set_ylim(bottom=0.58, top=1.02)
    gg.ax_joint.tick_params(axis='both', labelsize=12)
    gg.ax_joint.tick_params(
        axis='both', which='major', direction='inout', color='k',
        width=2, length=6, bottom=True, left=True, zorder=3
        )

    # set grid style
    gg.ax_joint.yaxis.grid(
        which="major", color='#d9d9d9', linestyle='--', linewidth=1
        )
    gg.ax_joint.xaxis.grid(
        which="major", color='#d9d9d9', linestyle='--', linewidth=1
        )
    gg.ax_joint.set_axisbelow(True)
    gg.ax_joint.add_artist(TickRedrawer())

    if c:
        # Plot mean and median
        _ = gg.ax_joint.axvline(
            x=df_stats['ani_mean'], ymin=0, ymax=1,
            color=vline_color, linewidth=2, linestyle='--',
            label='Mean'
            )
        _ = gg.ax_joint.axhline(
            y=df_stats['frag_mean'], xmin=0, xmax=1,
            color=vline_color, linewidth=2, linestyle='--',
            )
        _ = gg.ax_joint.axvline(
            x=df_stats['ani_median'], ymin=0, ymax=1,
            color=vline_color, linewidth=2, linestyle=':',
            label='Mean'
            )
        _ = gg.ax_joint.axhline(
            y=df_stats['frag_median'], xmin=0, xmax=1,
            color=vline_color, linewidth=2, linestyle=':',
            )

        # Build legend for mean and median
        gg.ax_joint.legend(
            loc='lower left',
            fontsize=12,
            markerscale=1.5,
            numpoints=1,
            frameon=False,
            ncol=2
            )

    # adjust layout, save, and close
    gg.fig.set_figwidth(7)
    gg.fig.set_figheight(5)
    gg.savefig(f'{outfile}_{species}.pdf')
    plt.close()


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--input_file',
        help='Please specify the input file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--output_file_prefix',
        help='Please specify the output file prefix!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-xmin', '--xaxis_minimum',
        help='OPTIONAL: Minimum value to plot on x-axis. (Default=95.0)',
        metavar='',
        type=float,
        default=95.0,
        required=False
        )
    parser.add_argument(
        '-xmax', '--xaxis_maximum',
        help='OPTIONAL: Maximum value to plot on x-axis. (Default=100.0)',
        metavar='',
        type=float,
        default=100.0,
        required=False
        )
    parser.add_argument(
        '-t', '--xaxis_step_size',
        help='OPTIONAL: X-axis ticks step increment. (Default=1.0)',
        metavar='',
        type=float,
        default=1.0,
        required=False
        )
    parser.add_argument(
        '-p', '--point_size',
        help='OPTIONAL: Size of the plotted points (Default=4.0)',
        metavar='',
        type=float,
        default=4.0,
        required=False
        )
    parser.add_argument(
        '-a', '--point_alpha',
        help='OPTIONAL: Alpha value of the plotted points (Default=0.10)',
        metavar='',
        type=float,
        default=0.10,
        required=False
        )
    parser.add_argument(
        '-c', '--add_cross_hairs',
        help='OPTIONAL: Input -c True mean/median cross hairs (Default=None).',
        metavar='',
        type=str,
        default=None,
        required=False
        )
    parser.add_argument(
        '-g', '--generate_GAM_trendline',
        help='OPTIONAL: Input -g True adds trendline with GAM (Default=None).',
        metavar='',
        type=str,
        default=None,
        required=False
        )
    parser.add_argument(
        '-r', '--random_subsample',
        help='OPTIONAL: Set > 1 to plot subsample of r genomes per species.',
        metavar='',
        type=int,
        default=1,
        required=False
        )
    parser.add_argument(
        '-e', '--repeat_subsamples',
        help='OPTIONAL: Repeat subsampling this many times (Default=100).',
        metavar='',
        type=int,
        default=100,
        required=False
        )
    parser.add_argument(
        '-s', '--single_species',
        help='OPTIONAL: Input -s True for single species plots (Default=None).',
        metavar='',
        type=str,
        default=None,
        required=False
        )
    parser.add_argument(
        '-l', '--all_species',
        help='OPTIONAL: Input -l True for all species one plot (Default=None).',
        metavar='',
        type=str,
        default=None,
        required=False
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\n\nRunning Script...\n')

    # define parameters
    infile = args['input_file']
    outfile = args['output_file_prefix']
    single_species = args['single_species']
    all_species = args['all_species']
    xmin = args['xaxis_minimum']
    xmax = args['xaxis_maximum']
    xstep = args['xaxis_step_size']
    p = args['point_size']
    a = args['point_alpha']
    z = None
    c = args['add_cross_hairs']
    g = args['generate_GAM_trendline']
    r = args['random_subsample']
    e = args['repeat_subsamples']

    # test params
    if not all_species and not single_species:
        print(
            '\nNo option specified. Please set -s, -a or both to True.\n\n'
            )
        sys.exit(1)

    # read in the data
    if r > 1:
        df, n = gather_subsampled_data(infile, xmin, xmax, r, e)
    else:
        df, n = gather_data(infile, xmin, xmax)

    # build the plot
    if all_species:
        if n > 10000: z = True
        else: z = None
        _ = fastANI_scatter_plot(
            df, n, 'All_species', outfile, xmin, xmax, xstep, p, a, z, g, c
            )
    if single_species:
        for species in df['species'].unique():
            print(f'\tPlotting {species} ...')
            dfx = df[df['species'] == species]
            n = len(dfx)
            if n > 10000: z = True
            else: z = None
            if n >= 20:
                _ = fastANI_scatter_plot(
                    dfx, n, species, outfile, xmin, xmax, xstep, p, a, z, g, c
                    )
                               
    print(f'\n\nComplete success space cadet!! Finished without errors.\n\n')

if __name__ == "__main__":
    main()
