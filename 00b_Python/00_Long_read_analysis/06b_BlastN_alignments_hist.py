 #!/usr/bin/env python

'''Returns a histogram of pID values from blastn.

pID = percent identity of the blast sequence alignment

For All vs All blast the long reads or contigs are set as the query and
as the reference database. This scripts needs tabular blast output with
this specific -outfmt flag order:

        -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart
        qend sstart send evalue bitscore qlen slen'

The aligned_frac is the (aligned query length) / (query length).

Default aligned_frac filter is >= 0.5 Wanted complete alignments for this.

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


def parse_blastn(infile):

    # Reads tabular blastn file, computes pid of alignments
    # returns array of floats

    alignment = []
    query = []
    fraction = []

    with open(infile, 'r') as file:
        for line in file:
            X = line.rstrip().split('\t')
            qname = X[0]
            rname = X[1]
            # check for self matches and discard
            if qname == rname: continue
            pid = float(X[2])
            alen = int(X[3]) # alignment length
            qlen = int(X[12]) # reference sequence length
            aligned_frac = alen / qlen

            if aligned_frac >= 0.1:
                alignment.append(alen)
                query.append(qlen)
                fraction.append(aligned_frac)

    return alignment, query, fraction


def build_hist_plot(alignment, query, fraction, bins, title, outfile):

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(7,10))

    fig.suptitle(title, fontsize=18, y=0.97)
    _ = ax1.hist(alignment, bins=bins, color='#c51b7d', alpha=0.6)
    ax1.set_ylabel('Read count', fontsize=12)
    ax1.set_xlabel('Alignment length', fontsize=12)

    _ = ax2.hist(query, bins=bins, color='#c51b7d', alpha=0.6)
    ax2.set_ylabel('Read count', fontsize=12)
    ax2.set_xlabel('Query read length', fontsize=12)

    _ = ax3.hist(fraction, bins=bins, color='#c51b7d', alpha=0.6)
    ax3.set_ylabel('Read count', fontsize=12)
    ax3.set_xlabel('Aligned fraction', fontsize=12)

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
        '-i', '--input_blastn_file',
        help='Please specify the input tabular BlastN file!',
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
    infile= args['input_blastn_file']
    outfile = args['output_file_name']
    title = args['plot_title']
    bins = args['number_of_bins']

    if not title:
        title = infile.split('.')[0]

    # parse the tabular BlastN file
    print('\nParsing the tabular BlastN file ...\n')
    alignment, query, fraction = parse_blastn(infile)

    # build some plots
    print('\nPlotting histogram ...\n')
    _ = build_hist_plot(alignment, query, fraction, bins, title, outfile)


    print(f'\n\nComplete success space cadet!! Finished without errors.\n\n')
    
if __name__ == "__main__":
    main()

