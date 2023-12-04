 #!/usr/bin/env python

'''Returns a histogram of percent sequence similarity from a PAF file.

Intended for reads mapped to a reference MAG or genome.
Generates two histograms of the percent sequence similarity distribution
of the mapped reads with x-axis 70-100 and 98-100. Generates histogram of
read length distribution.

Sequence similarity based on minimap2 manual:
https://lh3.github.io/minimap2/minimap2.html
https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity

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


def parse_paf(infile, min_pid, min_frac, min_qual):

    # Reads sam file, computes pid of alignments, returns array of floats

    data = {}
    total_match = 0
    dup_match = 0
    passing_match = 0

    with open(infile, 'r') as file:
        for line in file:
            total_match += 1
            X = line.split(sep='\t')
            qname = X[0] # query sequence name
            tname = X[5] # target sequence name
            if qname == tname: continue # remove any self matches
            qlen = int(X[1]) # length of the query sequence (read length)
            matches = int(X[9]) # number of matches across alignment length
            alen = int(X[10]) # length of the alignment
            pid = 100*(matches / alen)
            qual = int(X[11]) # mapping quality
            afrac = alen / qlen # aligned fractions
            # some quality and best match filtering
            # we aren't interesting in anything below 70% pid generally
            if aligned_frac >= min_frac and pid >= xmin and pid >= qscore:
                # count each query sequence only once.
                # For duplicate mathces keep the one with the highest qual
                if qname in data:
                    dup_match += 1
                    old_qual = data[qname][1]
                    if qual > old_qual:
                        data[qname] = [pid, qual, afrac, alen]
                else:
                    passing_match += 1
                    data[qname] = [pid, qual, afrac, alen]

    print(
        f'\n\t\tFile Name: {infile}\n'
        f'\t\tTotal matches in file: {total_match}\n'
        f'\t\tDuplicate matches: {dup_match}\n'
        f'\t\tTotal passing matches: {passing_match}\n'
        )

    return data


def build_pid_plots(pIDs, bins, title, outfile):

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
        '-i', '--input_paf_file',
        help='Please specify the PAF file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--output_file_prefix',
        help='Please specify the prefix to use for naming the output file!',
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
    parser.add_argument(
        '-min_frac', '--minimum_aligned_frac',
        help='(OPTIONAL) Specify the minimum aligned_frac (Default = 0.5).',
        metavar='',
        type=float,
        required=False,
        default=0.5
        )
    parser.add_argument(
        '-min_qual', '--minimum_quality',
        help='(OPTIONAL) Specify the minimum aligned_frac (Default = 5).',
        metavar='',
        type=int,
        required=False,
        default=5
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\n\nRunning Script...\n')

    # define parameters
    infile= args['input_paf_file']
    outpre = args['output_file_prefix']
    title = args['plot_title']
    bins = args['number_of_bins']
    xmin = args['minimum_value']
    min_frac = args['min_frac']
    min_qual = args['min_qual']

    if not title:
        title = infile.split('.')[0]

    # parse the sam file and compute percent sequence similarity (pID)
    print('\nParsing PAF file and computing percent sequence similarities ...\n')
    data = parse_paf(infile, min_pid, min_frac, min_qual)

    # build some plots
    print('\nPlotting histogram ...\n')
    pIDs = [v[0] for v in data.values()]
    _ = build_pid_plot(pIDs, bins, title, outpre)


    print(f'\n\nComplete success space cadet!! Finished without errors.\n\n')
    
if __name__ == "__main__":
    main()

