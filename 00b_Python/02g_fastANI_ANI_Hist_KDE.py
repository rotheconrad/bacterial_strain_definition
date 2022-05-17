 #!/usr/bin/env python

'''Two plots: Stepped histograms and KDEs for all species on each plot.

Each species is line showing the ANI distribution for genome pairs of
that species.

Color of the line correlates with number of genome pairs for the species.

    * Red (#e41a1c): >= 1,000,000 genome pairs
    * Blue (#377eb8): 100,00 - 999,999 genome pairs
    * Green (#4daf4a): 10,000 - 99,999 genome pairs
    * Purple (#984ea3): 1,000 - 9,999 genome pairs
    * Gray (#999999): 100 - 999 genome pairs

Outputs 3 files:
    - Count of genome pairs passing the filter for each species
    - histogram.pdf
    - kde.pdf

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

import argparse, sys
from collections import defaultdict
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import PercentFormatter
import seaborn as sns; sns.set(style="white", color_codes=True)
import numpy as np


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


def gather_data(ani_file, outfile, xmin, xmax):
    """Reads in the tsv ALLvsAll ANI """

    print("\nReading data.")

    # initialize dictionary to store {species: [ani]}
    ani_dict = defaultdict(list)
    # initialize dictionary to store {gname: tfrac}
    name_dict = defaultdict(dict)

    with open(ani_file, 'r') as f:
        for l in f:
            X = l.rstrip().split('\t')
            # get the genome pair names
            qry_genome = X[0]
            ref_genome = X[1]
            species = qry_genome.split('/')[1]
            # ani value
            ani = float(X[2])
            # filter xmin xmax
            if ani < xmin: continue
            elif ani > xmax: continue
            # total frac column. Keep larger genome as reference
            tfrac = int(X[4])
            # remove self matchs
            if qry_genome == ref_genome: continue
            # sort genome file names and combine
            names = [qry_genome, ref_genome]
            names.sort()
            gname = '-'.join(names)
            # record values from only one genome pair.
            # Keep larger genome as reference.
            if tfrac > name_dict[species].get(gname[0], 0):
                # add genome pair to dict
                name_dict[species][gname] = [tfrac, ani]

    # read through filtered genome pairs and store ANI by species
    for species, gpairs in name_dict.items():
        # store ani values for each species
        for gname, values in gpairs.items():
            ani = values[1]
            ani_dict[species].append(ani)

    print("\nWriting count of genome pairs passing the filter by species.")
    # write output file of species count
    with open(f'{outfile}_gpair_count.tsv', 'w') as ofile:
        for species, anis in ani_dict.items():
            ofile.write(f'{species}\t{len(anis)}\n')

    return ani_dict


def all_histogram_plot(ani_dict, outfile, xmin, xmax, xstep, hbins):

    # Set Colors and markers
    grid_color = '#d9d9d9'
    Red = '#e41a1c' # >= 1,000,000 genome pairs
    Blue = '#377eb8' # 100,00 - 999,999 genome pairs
    Green = '#4daf4a' # 10,000 - 99,999 genome pairs
    Purple = '#984ea3' #  1,000 - 9,999 genome pairs
    Gray = '#999999' # 100 - 999 genome pairs

    # Build the plot
    fig, (ax1, ax2) = plt.subplots(
        2, 1, gridspec_kw={'height_ratios': [4, 1]},
        )

    # plot the data
    count = 0 # species count
    failed = 0 # number not passing filter
    for species, ani in ani_dict.items():
        gpairs = len(ani)
        if gpairs >= 1000000:
            count += 1
            sns.histplot(
                x=ani, element="step", fill=False, stat='probability', 
                ax=ax1, bins=hbins, color=Red
                )
        elif gpairs >= 100000:
            count += 1
            sns.histplot(
                x=ani, element="step", fill=False, stat='probability', 
                ax=ax1, bins=hbins, color=Blue
                )
        elif gpairs >= 10000:
            count += 1
            sns.histplot(
                x=ani, element="step", fill=False, stat='probability', 
                ax=ax1, bins=hbins, color=Green
                )
        elif gpairs >= 1000:
            count += 1
            sns.histplot(
                x=ani, element="step", fill=False, stat='probability', 
                ax=ax1, bins=hbins, color=Purple
                )
        elif gpairs >= 100:
            count += 1
            sns.histplot(
                x=ani, element="step", fill=False, stat='probability', 
                ax=ax1, bins=hbins, color=Gray
                )
        else:
            failed += 1
            print(
                f'{species} not plotted. <100 genome pairs within x-range.'
                )

    print(f'\nSpecies failing filter: {failed}')
    print(f'Species passing filter: {count}')

    # plot title, labels, text
    ptitle = f'ANI Histogram for {count} species'
    fig.suptitle(ptitle, fontsize=18)

    ax1.set_xlabel(
        'Average nucleotide identity (%)',
        fontsize=12, y=-0.02
        )
    ax1.set_ylabel(
        'Probability',
        fontsize=12, x=-0.02
        )

    # set the axis parameters / style
    hstep = xstep/10
    ax1.set_xticks(np.arange(xmin, xmax+hstep, xstep))
    ax1.set_xlim(left=xmin-hstep, right=xmax+hstep)
    ax1.tick_params(axis='both', labelsize=12)
    ax1.tick_params(
        axis='both', which='major', direction='inout', color='k',
        width=2, length=6, bottom=True, left=True, zorder=3
        )

    # set grid style
    ax1.yaxis.grid(
        which="major", color=grid_color, linestyle='--', linewidth=1
        )
    ax1.xaxis.grid(
        which="major", color=grid_color, linestyle='--', linewidth=1
        )
    ax1.set_axisbelow(True)
    ax1.add_artist(TickRedrawer())

    # Build the Plot Legend
    reds = Line2D(
        [0],[0], color=Red, linestyle='-', lw=4,
        label='Species with >= 1,000,000 genomes pairs', 
        )
    blues = Line2D(
        [0],[0], color=Blue, linestyle='-', lw=4,
        label='Species with 100,00 - 999,999 genomes pairs', 
        )
    greens = Line2D(
        [0],[0], color=Green, linestyle='-', lw=4,
        label='Species with 10,000 - 99,999 genomes pairs', 
        )
    purples = Line2D(
        [0],[0], color=Purple, linestyle='-', lw=4,
        label='Species with 1,000 - 9,999 genomes pairs', 
        )
    grays = Line2D(
        [0],[0], color=Gray, linestyle='-', lw=4,
        label='Species with 100 - 999 genomes pairs', 
        )
    legend_elements = [
                        reds,
                        blues,
                        greens,
                        purples,
                        grays
                        ]
    ax2.legend(
        handles=legend_elements,
        loc='center',
        fontsize=12,
        fancybox=True,
        framealpha=0.0,
        frameon=False,
        ncol=1
        )
    ax2.axis('off')

    # adjust layout, save, and close
    fig.set_figwidth(7)
    fig.set_figheight(5)
    fig.tight_layout()
    fig.savefig(f'{outfile}_hist_All.pdf')
    plt.close()

    return True


def all_KDE_plot(ani_dict, outfile, xmin, xmax, xstep, hbins):
    
    # Set Colors and markers
    grid_color = '#d9d9d9'
    Red = '#e41a1c' # >= 1,000,000 genome pairs
    Blue = '#377eb8' # 100,00 - 999,999 genome pairs
    Green = '#4daf4a' # 10,000 - 99,999 genome pairs
    Purple = '#984ea3' #  1,000 - 9,999 genome pairs
    Gray = '#999999' # 100 - 999 genome pairs

    # Build the plot
    fig, (ax1, ax2) = plt.subplots(
        2, 1, gridspec_kw={'height_ratios': [4, 1]},
        )

    # plot the data
    count = 0 # species count
    failed = 0 # number not passing filter
    for species, ani in ani_dict.items():
        gpairs = len(ani)
        if gpairs >= 1000000:
            count += 1
            sns.kdeplot(
                x=ani, ax=ax1, color=Red, common_norm=True
                )
        elif gpairs >= 100000:
            count += 1
            sns.kdeplot(
                x=ani, ax=ax1, color=Blue, common_norm=True
                )
        elif gpairs == 16125: #>= 10000:
            count += 1
            sns.kdeplot(
                x=ani, ax=ax1, color=Green, common_norm=True
                )
        elif gpairs >= 1000:
            count += 1
            sns.kdeplot(
                x=ani, ax=ax1, color=Purple, common_norm=True
                )
        elif gpairs >= 100:
            count += 1
            sns.kdeplot(
                x=ani, ax=ax1, color=Gray, common_norm=True
                )
        else:
            failed += 1
            print(
                f'{species} not plotted. <100 genome pairs within x-range.'
                )

    print(f'\nSpecies failing filter: {failed}')
    print(f'Species passing filter: {count}')

    # Set KDE to probability instead of density
    binwidth = (xmax - xmin) / hbins
    ax1.yaxis.set_major_formatter(PercentFormatter(1 / binwidth))

    # plot title, labels, text
    ptitle = f'ANI KDE for {count} species'
    fig.suptitle(ptitle, fontsize=18)

    ax1.set_xlabel(
        'Average nucleotide identity (%)',
        fontsize=12, y=-0.02
        )
    ax1.set_ylabel(
        'Probability',
        fontsize=12, x=-0.02
        )

    # set the axis parameters / style
    hstep = xstep/10
    ax1.set_xticks(np.arange(xmin, xmax+hstep, xstep))
    ax1.set_xlim(left=xmin-hstep, right=xmax+hstep)
    ax1.tick_params(axis='both', labelsize=12)
    ax1.tick_params(
        axis='both', which='major', direction='inout', color='k',
        width=2, length=6, bottom=True, left=True, zorder=3
        )

    # set grid style
    ax1.yaxis.grid(
        which="major", color=grid_color, linestyle='--', linewidth=1
        )
    ax1.xaxis.grid(
        which="major", color=grid_color, linestyle='--', linewidth=1
        )
    ax1.set_axisbelow(True)
    ax1.add_artist(TickRedrawer())

    # Build the Plot Legend
    reds = Line2D(
        [0],[0], color=Red, linestyle='-', lw=4,
        label='Species with >= 1,000,000 genomes pairs', 
        )
    blues = Line2D(
        [0],[0], color=Blue, linestyle='-', lw=4,
        label='Species with 100,00 - 999,999 genomes pairs', 
        )
    greens = Line2D(
        [0],[0], color=Green, linestyle='-', lw=4,
        label='Species with 10,000 - 99,999 genomes pairs', 
        )
    purples = Line2D(
        [0],[0], color=Purple, linestyle='-', lw=4,
        label='Species with 1,000 - 9,999 genomes pairs', 
        )
    grays = Line2D(
        [0],[0], color=Gray, linestyle='-', lw=4,
        label='Species with 100 - 999 genomes pairs', 
        )

    legend_elements = [
                        reds,
                        blues,
                        greens,
                        purples,
                        grays
                        ]

    ax2.legend(
        handles=legend_elements,
        loc='center',
        fontsize=12,
        fancybox=True,
        framealpha=0.0,
        frameon=False,
        ncol=1
        )
    ax2.axis('off')

    # adjust layout, save, and close
    fig.set_figwidth(7)
    fig.set_figheight(5)
    fig.tight_layout()
    fig.savefig(f'{outfile}_KDE_All.pdf')
    plt.close()

    return True


def single_histogram_plot(species, ani, outfile, xmin, xmax, xstep, hbins):

    # Set Colors and markers
    grid_color = '#d9d9d9'

    # Build the plot
    fig, (ax1, ax2) = plt.subplots(
        2, 1, gridspec_kw={'height_ratios': [9, 1]},
        )

    # plot the data
    gpairs = len(ani)

    sns.histplot(
        x=ani, element="step", fill=False, stat='probability', 
        ax=ax1, bins=hbins, color='k'
        )

    # plot title, labels, text
    ptitle = f'ANI Histogram for {species}'
    fig.suptitle(ptitle, fontsize=18)

    ax1.set_xlabel(
        'Average nucleotide identity (%)',
        fontsize=12, y=-0.02
        )
    ax1.set_ylabel(
        'Probability',
        fontsize=12, x=-0.02
        )
    ax2.set_title(
        f'Genome pairs [{xmin}, {xmax}]: {gpairs}',
        fontsize=18, #y=-0.5
        )

    # set the axis parameters / style
    hstep = xstep/10
    ax1.set_xticks(np.arange(xmin, xmax+hstep, xstep))
    ax1.set_xlim(left=xmin-hstep, right=xmax+hstep)
    ax1.tick_params(axis='both', labelsize=12)
    ax1.tick_params(
        axis='both', which='major', direction='inout', color='k',
        width=2, length=6, bottom=True, left=True, zorder=3
        )

    # set grid style
    ax1.yaxis.grid(
        which="major", color=grid_color, linestyle='--', linewidth=1
        )
    ax1.xaxis.grid(
        which="major", color=grid_color, linestyle='--', linewidth=1
        )
    ax1.set_axisbelow(True)
    ax1.add_artist(TickRedrawer())

    # adjust layout, save, and close
    ax2.axis('off')
    fig.set_figwidth(7)
    fig.set_figheight(5)
    fig.tight_layout()
    fig.savefig(f'{outfile}_hist_{species}.pdf')
    plt.close()

    return True


def single_KDE_plot(species, ani, outfile, xmin, xmax, xstep, hbins):
    
    # Set Colors and markers
    grid_color = '#d9d9d9'

    # Build the plot
    fig, (ax1, ax2) = plt.subplots(
        2, 1, gridspec_kw={'height_ratios': [9, 1]},
        )

    # plot the data
    gpairs = len(ani)
    sns.kdeplot(x=ani, ax=ax1, color='k')

    # Set KDE to probability instead of density
    binwidth = (xmax - xmin) / hbins
    ax1.yaxis.set_major_formatter(PercentFormatter(1 / binwidth))

    # plot title, labels, text
    ptitle = f'ANI KDE for {species}'
    fig.suptitle(ptitle, fontsize=18)

    ax1.set_xlabel(
        'Average nucleotide identity (%)',
        fontsize=12, y=-0.02
        )
    ax1.set_ylabel(
        'Probability',
        fontsize=12, x=-0.02
        )
    ax2.set_title(
        f'Genome pairs between [{xmin}, {xmax}]: {gpairs}',
        fontsize=18, #y=-0.5
        )
    # set the axis parameters / style
    hstep = xstep/10
    ax1.set_xticks(np.arange(xmin, xmax+hstep, xstep))
    ax1.set_xlim(left=xmin-hstep, right=xmax+hstep)
    ax1.tick_params(axis='both', labelsize=12)
    ax1.tick_params(
        axis='both', which='major', direction='inout', color='k',
        width=2, length=6, bottom=True, left=True, zorder=3
        )
    # set grid style
    ax1.yaxis.grid(
        which="major", color=grid_color, linestyle='--', linewidth=1
        )
    ax1.xaxis.grid(
        which="major", color=grid_color, linestyle='--', linewidth=1
        )
    ax1.set_axisbelow(True)
    ax1.add_artist(TickRedrawer())

    # adjust layout, save, and close
    ax2.axis('off')
    fig.set_figwidth(7)
    fig.set_figheight(5)
    fig.tight_layout()
    fig.savefig(f'{outfile}_KDE_{species}.pdf')
    plt.close()

    return True


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
        help='Please specify the output file prefix?',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-hist', '--plot_histogram',
        help='OPTIONAL: Input -hist True to plot histograms (Default=None).',
        metavar='',
        type=str,
        default=None,
        required=False
        )
    parser.add_argument(
        '-kde', '--plot_kde',
        help='OPTIONAL: Input -kde True to plot KDEs (Default=None).',
        metavar='',
        type=str,
        default=None,
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
        '-a', '--all_species',
        help='OPTIONAL: Input -a True for all species one plot (Default=None).',
        metavar='',
        type=str,
        default=None,
        required=False
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
        '-b', '--hist_bins',
        help='OPTIONAL: Number of bins for histogram. (Default=30)',
        metavar='',
        type=int,
        default=30,
        required=False
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\n\nRunning Script...\n')

    # define parameters
    infile = args['input_file']
    outfile = args['output_file_prefix']
    hist = args['plot_histogram']
    kde = args['plot_kde']
    single_species = args['single_species']
    all_species = args['all_species']
    xmin = args['xaxis_minimum']
    xmax = args['xaxis_maximum']
    xstep = args['xaxis_step_size']
    hbins = args['hist_bins']

    # make the plots
    if not hist and not kde:
        print(
            '\nNo option specified. Please set -hist, -kde or both to True.\n\n'
            )
        sys.exit(1)

    if not all_species and not single_species:
        print(
            '\nNo option specified. Please set -s, -a or both to True.\n\n'
            )
        sys.exit(1)

    # read in the data
    ani_dict = gather_data(infile, outfile, xmin, xmax)

    if all_species:
        if hist:
            _ = all_histogram_plot(
                                ani_dict, outfile, xmin, xmax, xstep, hbins
                                )
        if kde:
            _ = all_KDE_plot(ani_dict, outfile, xmin, xmax, xstep, hbins)

    if single_species:
        if hist:
            for species, ani in ani_dict.items():
                print(f'Plotting histogram for {species}')
                _ = single_histogram_plot(
                                species, ani, outfile, xmin, xmax, xstep, hbins
                                )
        if kde:
            for species, ani in ani_dict.items():
                print(f'Plotting KDE for {species}')
                _ = single_KDE_plot(
                                species, ani, outfile, xmin, xmax, xstep, hbins
                                )

    print(f'\n\nComplete success space cadet!! Finished without errors.\n\n')

if __name__ == "__main__":
    main()
