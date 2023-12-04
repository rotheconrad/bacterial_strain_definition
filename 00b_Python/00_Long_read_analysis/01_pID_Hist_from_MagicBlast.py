 #!/usr/bin/env python

'''Returns a histogram of percent sequence similarity from magicblast.

Intended for reads mapped to a reference MAG or genome.
Intended for tabular magic blast output that has already been filtered.
Generates a histogram of the percent sequence similarity distribution
of the mapped reads.

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: September 2022
License :: GNU GPLv3
Copyright 2022 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse
import matplotlib.pyplot as plt


def parse_magicblast(infile, xmin):

    # Reads sam file, computes pid of alignments, returns array of floats

    pIDs = []

    with open(infile, 'r') as file:
        for line in file:
            if line.startswith('#'): continue
            X = line.rstrip().split('\t')
            pid = float(X[2])
            if pid >= xmin: pIDs.append(pid)

    return pIDs


def build_hist_plot(pIDs, bins, title, outfile):

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7,5))

    fig.suptitle(title, fontsize=18, y=0.95)
    _ = ax1.hist(pIDs, bins=bins, color='#c51b7d', alpha=0.6)
    ax1.set_ylabel('Alignment count', fontsize=12)
    ax1.text(
            0.05, 0.90, f'n = {len(pIDs)}',
            fontsize=8, transform=ax1.transAxes
            )

    pIDs98 = [i for i in pIDs if i >= 98]
    _ = ax2.hist(pIDs98, bins=bins, color='#c51b7d', alpha=0.6)
    ax2.set_ylabel('Alignment count', fontsize=12)
    ax2.set_xlabel('Sequence alignment (%)', fontsize=12)
    ax2.text(
            0.05, 0.90, f'n = {len(pIDs98)}',
            fontsize=8, transform=ax2.transAxes
            )

    for ax in [ax1, ax2]:
        ax.minorticks_on()
        ax.tick_params(
            which='minor', axis='both', left=False, bottom=False
            )
        ax.tick_params(
                    which='major', axis='both',
                    left=True, bottom=True,
                    size=6, width=2, tickdir='inout',
                    labelsize=12, zorder=10
                    )
        ax.yaxis.grid(
            which="major", color='#bdbdbd', linestyle='--',
            linewidth=1, alpha=0.4, zorder=1
            )
        ax.set_axisbelow(True)
        for spine in ax.spines.values(): spine.set_linewidth(2)

    #fig.set_tight_layout(True)
    plt.subplots_adjust(left=0.15, right=0.98)
    plt.savefig(outfile)
    plt.close()

    return True


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--input_magicblast_file',
        help='Please specify the filtered tabular magicblast file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--output_file_name',
        help='Please specify the output file name (use .pdf)!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-t', '--plot_title',
        help='(OPTIONAL) Please specify the plot title!',
        metavar='',
        type=str,
        required=False
        )
    parser.add_argument(
        '-b', '--number_of_bins',
        help='(OPTIONAL) Specify the number of bins (Default = 30).',
        metavar='',
        type=str,
        required=False,
        default=30
        )
    parser.add_argument(
        '-xmin', '--minimum_value',
        help='(OPTIONAL) Specify the miniman value to include (Default = 70).',
        metavar='',
        type=str,
        required=False,
        default=70
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\n\nRunning Script...\n')

    # define parameters
    infile= args['input_magicblast_file']
    outfile = args['output_file_name']
    title = args['plot_title']
    bins = args['number_of_bins']
    xmin = args['minimum_value']

    if not title:
        title = infile.split('.')[0]

    # parse the magicblast file and grab the percent sequence similarity (pID)
    print('\nParsing MagicBlast file to grab percent sequence similarities ...\n')
    pIDs = parse_magicblast(infile, xmin)

    # build some plots
    print('\nPlotting histogram ...\n')
    _ = build_hist_plot(pIDs, bins, title, outfile)


    print(f'\n\nComplete success space cadet!! Finished without errors.\n\n')
    
if __name__ == "__main__":
    main()

