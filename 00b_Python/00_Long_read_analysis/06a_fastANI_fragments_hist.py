 #!/usr/bin/env python

'''Returns a histogram of ANI values from fastANI.

Filters by shared fraction of the sequence pair.

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


def parse_fastani(infile):

    # Reads fastani file, returns array of floats

    shared = []
    total = []
    fraction = []

    with open(infile, 'r') as file:
        for line in file:
            X = line.rstrip().split('\t')
            qry = X[0]
            ref = X[1]
            # check for self matches and discard
            if qry == ref: continue
            ani = float(X[2])
            shared_frags = int(X[3])
            total_frags = int(X[4])
            shared_fraction = shared_frags / total_frags

            shared.append(shared_frags)
            total.append(total_frags)
            fraction.append(shared_fraction)
    
    return shared, total, fraction


def build_hist_plot(shared, total, fraction, bins, title, outfile):

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(7,10))

    fig.suptitle(title, fontsize=18, y=0.97)
    _ = ax1.hist(shared, bins=bins, color='#c51b7d', alpha=0.6)
    ax1.set_ylabel('Read count', fontsize=12)
    ax1.set_xlabel('Shared fragments', fontsize=12)

    _ = ax2.hist(total, bins=bins, color='#c51b7d', alpha=0.6)
    ax2.set_ylabel('Read count', fontsize=12)
    ax2.set_xlabel('Total fragments', fontsize=12)

    _ = ax3.hist(fraction, bins=bins, color='#c51b7d', alpha=0.6)
    ax3.set_ylabel('Read count', fontsize=12)
    ax3.set_xlabel('Shared fraction', fontsize=12)

    for ax in [ax1, ax2, ax3]:
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
        '-i', '--input_fastANI_file',
        help='Please specify the input all vs all fastANI file!',
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
        type=int,
        required=False,
        default=30
        )

    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\n\nRunning Script...\n')

    # define parameters
    infile= args['input_fastANI_file']
    outfile = args['output_file_name']
    title = args['plot_title']
    bins = args['number_of_bins']

    if not title:
        title = infile.split('.')[0]

    # parse the fastANI file
    print('\nParsing fastANI file ...\n')
    shared, total, fraction = parse_fastani(infile)

    # build some plots
    print('\nPlotting histogram ...\n')
    _ = build_hist_plot(shared, total, fraction, bins, title, outfile)


    print(f'\n\nComplete success space cadet!! Finished without errors.\n\n')
    
if __name__ == "__main__":
    main()

