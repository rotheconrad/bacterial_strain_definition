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


def parse_blastn(infile, min_pid, min_frac):

    # Reads tabular blastn file, computes pid of alignments
    # returns array of floats

    pid_dict = {}
    total_match = 0
    self_match = 0
    reciprocal_match = 0
    passing_match = 0

    with open(infile, 'r') as file:
        for line in file:
            total_match += 1
            X = line.rstrip().split('\t')
            qname = X[0]
            rname = X[1]
            # check for self matches and discard
            if qname == rname:
                self_match += 1
                continue
            pid = float(X[2])
            bitscore = float(X[11])
            alen = int(X[3]) # alignment length
            qlen = int(X[12]) # reference sequence length
            aligned_frac = alen / qlen

            if aligned_frac >= min_frac and pid >= min_pid:
                # check for reciprocal matches and discard lower score
                sorted_pair = '-'.join(sorted([qname, rname]))
                if sorted_pair in pid_dict:
                    reciprocal_match += 1
                    old_bitscore = pid_dict[sorted_pair][1]
                    if bitscore > old_bitscore:
                        pid_dict[sorted_pair] = [pid, bitscore]
                else:
                    passing_match += 1
                    pid_dict[sorted_pair] = [pid, bitscore]
    
    pIDs = [v[0] for v in pid_dict.values()]

    print(
        f'\n\t\tFile Name: {infile}\n'
        f'\t\tTotal matches in file: {total_match}\n'
        f'\t\tSelf matches: {self_match}\n'
        f'\t\tReciprocal matches: {reciprocal_match}\n'
        f'\t\tTotal passing matches: {passing_match}\n'
        )

    return pIDs


def build_hist_plot(pIDs, bins, title, outpre):

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7,5))

    fig.suptitle(title, fontsize=18, y=0.97)
    _ = ax1.hist(pIDs, bins=bins, color='#c51b7d', alpha=0.6)
    ax1.set_ylabel('Read count', fontsize=12)
    ax1.text(
            0.05, 0.90, f'n = {len(pIDs)}',
            fontsize=8, transform=ax1.transAxes
            )

    pIDs98 = [i for i in pIDs if i >= 98]
    _ = ax2.hist(pIDs98, bins=bins, color='#c51b7d', alpha=0.6)
    ax2.set_ylabel('Read count', fontsize=12)
    ax2.set_xlabel('Sequence identity (%)', fontsize=12)
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
    plt.savefig(f'{outpre}.pdf')
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
    parser.add_argument(
        '-min_pid', '--minimum_pid_value',
        help='(OPTIONAL) Specify the minimum pID value (Default = 70).',
        metavar='',
        type=int,
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
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\n\nRunning Script...\n')

    # define parameters
    infile= args['input_blastn_file']
    outfile = args['output_file_name']
    title = args['plot_title']
    bins = args['number_of_bins']
    min_pid = args['minimum_pid_value']
    min_frac = args['minimum_aligned_frac']

    if not title:
        title = infile.split('.')[0]

    # parse the tabular BlastN file
    print('\nParsing the tabular BlastN file ...\n')
    pIDs = parse_blastn(infile, min_pid, min_frac)

    # build some plots
    print('\nPlotting histogram ...\n')
    _ = build_hist_plot(pIDs, bins, title, outfile)


    print(f'\n\nComplete success space cadet!! Finished without errors.\n\n')
    
if __name__ == "__main__":
    main()

