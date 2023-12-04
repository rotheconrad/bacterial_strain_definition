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


def parse_fastani(infile, min_ani, min_frac):

    # Reads fastani file, computes pid of alignments, returns array of floats

    ani_dict = {}
    total_match = 0
    self_match = 0
    reciprocal_match = 0
    passing_match = 0

    with open(infile, 'r') as file:
        for line in file:
            total_match += 1
            X = line.rstrip().split('\t')
            qry = X[0]
            ref = X[1]
            # check for self matches and discard
            if qry == ref:
                self_match += 1
                continue
            ani = float(X[2])
            shared_frags = int(X[3])
            total_frags = int(X[4])
            shared_fraction = shared_frags / total_frags

            if shared_fraction >= min_frac and ani >= min_ani:
                # check for reciprocal matches and discard lower score
                sorted_pair = '-'.join(sorted([qry, ref]))
                if sorted_pair in ani_dict:
                    reciprocal_match += 1
                    old_ani = ani_dict[sorted_pair]
                    ani_dict[sorted_pair] = max([ani, old_ani])
                else:
                    passing_match += 1
                    ani_dict[sorted_pair] = ani
    
    ANIs = [i for i in ani_dict.values()]

    print(
        f'\n\t\tFile Name: {infile}\n'
        f'\t\tTotal matches in file: {total_match}\n'
        f'\t\tSelf matches: {self_match}\n'
        f'\t\tReciprocal matches: {reciprocal_match}\n'
        f'\t\tTotal passing matches: {passing_match}\n'
        )

    return ANIs


def build_hist_plot(ANIs, bins, title, outfile):

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7,5))

    fig.suptitle(title, fontsize=18, y=0.97)
    _ = ax1.hist(ANIs, bins=bins, color='#c51b7d', alpha=0.6)
    ax1.set_ylabel('Read count', fontsize=12)
    ax1.text(
            0.05, 0.90, f'n = {len(ANIs)}',
            fontsize=8, transform=ax1.transAxes
            )

    ANIs98 = [i for i in ANIs if i >= 98]
    _ = ax2.hist(ANIs98, bins=bins, color='#c51b7d', alpha=0.6)
    ax2.set_ylabel('Read count', fontsize=12)
    ax2.set_xlabel('Average nucleotide identity (%)', fontsize=12)
    ax2.text(
            0.05, 0.90, f'n = {len(ANIs98)}',
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
    parser.add_argument(
        '-min_ani', '--minimum_ani_value',
        help='(OPTIONAL) Specify the minimum ani value (Default = 70).',
        metavar='',
        type=int,
        required=False,
        default=70
        )
    parser.add_argument(
        '-min_frac', '--minimum_frac_value',
        help='(OPTIONAL) Specify the minimum frac value (Default = 0.5).',
        metavar='',
        type=float,
        required=False,
        default=0.5
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\n\nRunning Script...\n')

    # define parameters
    infile= args['input_fastANI_file']
    outfile = args['output_file_name']
    title = args['plot_title']
    bins = args['number_of_bins']
    min_ani = args['minimum_ani_value']
    min_frac = args['minimum_frac_value']

    if not title:
        title = infile.split('.')[0]

    # parse the fastANI file
    print('\nParsing fastANI file ...\n')
    ANIs = parse_fastani(infile, min_ani, min_frac)

    # build some plots
    print('\nPlotting histogram ...\n')
    _ = build_hist_plot(ANIs, bins, title, outfile)


    print(f'\n\nComplete success space cadet!! Finished without errors.\n\n')
    
if __name__ == "__main__":
    main()

