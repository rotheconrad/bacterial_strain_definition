#!/usr/bin/env python

'''Plots a histogram of all Phred scores in a fastq file.

Inputs:
    1) a fastq file to evaluate
    2) the Phred.tsv file linking ascii characters to phred values.

Outputs:
    1) a PDF file of histogram plot.

!!! Requires Phred.tsv file - located in github repo !!!

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: Oct 2022
License :: GNU GPLv3
Copyright 2022 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse
import matplotlib.pyplot as plt

def make_phred_table(qfile):
    d = {}
    with open(qfile, 'r') as p:
        for l in p:
            X = l.rstrip().split('\t')
            Q = int(X[0])
            C = X[1]
            d[C] = Q

    return d


def read_fastq(fq):
    linecount = 0
    while fq:
        name = fq.readline().rstrip().split(' ')[0][1:]
        if not name: break
        seq = fq.readline().rstrip()
        blank = fq.readline().rstrip()
        qal = fq.readline().rstrip()
        linecount += 4
        if linecount % 4 == 0: yield (name, seq, qal)


def parse_fastq_file(infile, phred):
    # reads fastq files. Convert ascii quality score to phred value.
    # return two array of read length and phred values

    readlen = []
    qscores = []

    with open(infile, 'r') as file:
        for name, seq, qal in read_fastq(file):
            readlen.append(len(seq))
            for q in qal:
                qscores.append(phred[q])

    return readlen, qscores


def build_hist_plot(readlen, qscores, bins, title, outfile):

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7,5))

    _ = ax1.hist(readlen, bins=bins, color='#c51b7d', alpha=0.6)
    ax1.set_title('Distribution of Read Lengths')
    ax1.set_xlabel('Read length (bp)')
    ax1.set_ylabel('Count')

    _ = ax2.hist(qscores, bins=bins, color='#c51b7d', alpha=0.6)
    ax2.set_title('Distribution of Phred Scores')
    ax2.set_xlabel('Phred score')
    ax2.set_ylabel('Count')

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
        '-i', '--input_fastq_file',
        help='Please specify the input fastq files!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-p', '--Phred_tsv_file',
        help='Please specify the Phred.tsv file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--out_file_name',
        help='Please specify the output file name!',
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
    infile= args['input_fastq_file']
    qfile = args['Phred_tsv_file']
    outfile = args['out_file_name']
    title = args['plot_title']
    bins = args['number_of_bins']

    if not title:
        title = infile.split('.')[0]

    # Retrieve Phred Score ASCII characters
    print('\nParsing the Phred.tsv file ...\n')
    phred = make_phred_table(args['Phred_tsv_file'])

    # parse the tabular BlastN file
    print('\nRetrieving read length & quality scores from fastq file ...\n')
    readlen, qscores = parse_fastq_file(infile, phred)

    # build some plots
    print('\nPlotting histogram ...\n')
    _ = build_hist_plot(readlen, qscores, bins, title, outfile)

    print(f'\n\nComplete success space cadet!! Finished without errors.\n\n')

if __name__ == "__main__":
    main()