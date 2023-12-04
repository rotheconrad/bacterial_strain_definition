 #!/usr/bin/env python

''' Compare all vs all results from fastANI to BlastN

Reads and filters equivalent all vs all fastANI and BlastN tabular 
output files. Calculates the differences between the scores of pairwise
comparisons.

Requires a read pair to file to key for fastANI and Blast ANI to BlastN.
BlastANI and FastANI were built for genomes, so
each sequence to compare needs its own fasta file and the output file
contains file names for query and reference whereas BlastN accepts and 
compares multiple sequences in a fasta file and the default read names
are used for the query and reference. We need them to match. I created
a two column text file separated with a single white space. The file name
(and the path if it is in the fastANI output file) in column 1 and the
corresponding read name (this should be the fasta defline) in column 2.
This script expects column two to have the '>' character of fasta files.

ex:
directory/path/read_014437.fasta >m64078_201009_211857

Outputs tsv file with ANI, pID and difference for each
pairwise comparison found in both results.

Outputs results from fastANI only.

Outputs results from BlastN only.

Output header key:
    - Read pair - name of the matched read pair
    - x diff: ANI - pID
    - y diff: shared fraction - aligned fraction
    - ANI: Score from fastANI
    - pID: Score from BlastN
    - sfrac: shared fraction
    - afrac: aligned fraction

    * pID is the percent sequence alignment from BlastN
    * shared fraction = shared fragments / total fragments
    * aligned fraction = alignment length / query length

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


def file_to_read_key(key_file):
    # creates a dict from the two column input file
    # {file-name: read-name}
    # returns dict

    fast_key = {}

    with open(key_file, 'r') as file:
        for line in file:
            X = line.rstrip().split(' ')
            fname = X[0]
            rname = X[1][1:]
            fast_key[fname] = rname

    return fast_key


def parse_blastani(infile, fast_key):

    blastani_dict = {}

    with open(infile, 'r') as file:
        header = file.readline()
        for line in file:
            X = line.rstrip().split('\t')
            qry = fast_key[X[0].split(':')[0] + ".fasta"]
            ref = fast_key[X[0].split(':')[1] + ".fasta"]
            ani = float(X[1])
            frac = float(X[2])

            blastani_dict[f"{qry}-{ref}"] = [ani, frac]

    return blastani_dict


def parse_fastani(infile, min_ani, min_frac, fast_key):

    # Reads fastani file, filters results, returns dict

    # {read-pair: [ani, shared_frac]}
    fast_dict = {}

    with open(infile, 'r') as file:
        for line in file:
            X = line.rstrip().split('\t')
            qry = fast_key[X[0]]
            ref = fast_key[X[1]]
            # check for self matches and discard
            if qry == ref: continue
            ani = float(X[2])
            shared_frags = int(X[3])
            total_frags = int(X[4])
            shared_frac = shared_frags / total_frags

            if shared_frac >= min_frac and ani >= min_ani:
                # check for reciprocal matches and discard lower score
                sorted_pair = '-'.join(sorted([qry, ref]))
                if sorted_pair in fast_dict:
                    old_ani = fast_dict[sorted_pair][0]
                    if ani > old_ani:
                        fast_dict[sorted_pair] = [ani, shared_frac]
                else:
                    fast_dict[sorted_pair] = [ani, shared_frac]

    return fast_dict


def parse_blastn(infile, min_pid, min_frac):

    # Reads tabular blastn file, filters results, returns dict

    # {read-pair: [pid, aligned_frac]}
    blast_dict = {}

    with open(infile, 'r') as file:
        for line in file:
            X = line.rstrip().split('\t')
            qname = X[0]
            rname = X[1]
            # check for self matches and discard
            if qname == rname: continue
            pid = float(X[2])
            bitscore = float(X[11])
            alen = int(X[3]) # alignment length
            qlen = int(X[12]) # reference sequence length
            aligned_frac = alen / qlen

            if aligned_frac >= min_frac and pid >= min_pid:
                # check for reciprocal matches and discard lower score
                sorted_pair = '-'.join(sorted([qname, rname]))
                if sorted_pair in blast_dict:
                    old_bitscore = blast_dict[sorted_pair][2]
                    if bitscore > old_bitscore:
                        blast_dict[sorted_pair] = [pid, aligned_frac, bitscore]
                else:
                    blast_dict[sorted_pair] = [pid, aligned_frac, bitscore]

    return blast_dict


def compute_difference(blastANI, fastANI, BlastN):

    # read blastANI matches and compute difference to fastANI and BlastN

    A = [] # blastANI vs fastANI, ANI
    B = [] # blastANI vs BlastN, ANI
    C = [] # blastANI vs fastANI, fragments
    D = [] # blastANI vs BlastN, fragments

    for k, v in blastANI.items():
        # blast ANI
        bani = v[0]
        bfrac = v[1]
        # fast ANI

        try:
            fani = fastANI[k][0]
            A.append(bani - fani)
        except: pass
        try:
            ffrac = fastANI[k][1]
            B.append(bani - pid)
        except: pass
        # blast N
        try:
            pid = BlastN[k][0]
            C.append(bfrac - ffrac)
        except: pass
        try:
            aln = BlastN[k][1]
            D.append(bfrac - aln)
        except: pass

    data = {"A": A, "B": B, "C": C, "D": D}

    return data


def build_plots(data, bins, title, outfile):

    # A blastANI vs fastANI, ANI
    # B blastANI vs BlastN, ANI
    # C blastANI vs fastANI, fragments
    # D blastANI vs BlastN, fragments

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(7,7))

    fig.suptitle(title, fontsize=18, y=0.97)

    _ = ax1.hist(data["A"], bins=bins, color='#c51b7d', alpha=0.6)
    ax1.set_xlabel('ANI diff (blastANI - fastANI)', fontsize=10)
    ax1.set_ylabel('Count', fontsize=10)

    _ = ax2.hist(data["B"], bins=bins, color='#c51b7d', alpha=0.6)
    ax2.set_xlabel('ANI diff (blastANI - blastN)', fontsize=10)
    ax2.set_ylabel('Count', fontsize=10)

    _ = ax3.hist(data["C"], bins=bins, color='#c51b7d', alpha=0.6)
    ax3.set_xlabel('Frag diff (blastANI - fastANI)', fontsize=10)
    ax3.set_ylabel('Count', fontsize=10)

    _ = ax4.hist(data["D"], bins=bins, color='#c51b7d', alpha=0.6)
    ax4.set_xlabel('Frag diff (blastANI - blastN)', fontsize=10)
    ax4.set_ylabel('Count', fontsize=10)

    for ax in [ax1, ax2, ax3, ax4]:
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
        '-ba', '--input_blastani_file',
        help='Please specify the input Blast ANI file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-fa', '--input_fastani_file',
        help='Please specify the input fastANI file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-k', '--fastani_file_key',
        help='Please specify file to read name key for fastANI!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-bn', '--input_blastn_file',
        help='Please specify the input tabular BlastN file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--output_file_name',
        help='Please specify the output file name (Use .pdf)!',
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
    ANIblastin = args['input_blastani_file']
    fastin = args['input_fastani_file']
    key_file = args['fastani_file_key']
    blastin = args['input_blastn_file']
    outfile = args['output_file_name']
    title = args['plot_title']
    bins = args['number_of_bins']
    min_ani = args['minimum_ani_value']
    min_pid = args['minimum_pid_value']
    min_frac = args['minimum_aligned_frac']

    if not title:
        title = "Blast ANI vs fastANI or BlastN"

    # get the read name to file key for fastANI
    print('\nParsing file name to read name key ...\n')
    fast_key = file_to_read_key(key_file)

    # parse the blastANI file
    print('\nParsing Blast ANI file ...\n')
    blastANI = parse_blastani(ANIblastin, fast_key)

    # parse the fastANI file - returns dict {read_pair: [ANI, shared_fraction]}
    print('\nParsing fastANI file ...\n')
    fastANI = parse_fastani(fastin, min_ani, min_frac, fast_key)

    # parse the BlastN file - returns dict {read_pair: [pID, aligned_fraction]}
    print('\nParsing the tabular BlastN file ...\n')
    BlastN = parse_blastn(blastin, min_pid, min_frac)

    # compare the two
    data = compute_difference(blastANI, fastANI, BlastN)

    # build some plots
    print('\nBuilding plots ...')
    _ = build_plots(data, bins, title, outfile)

    print(f'\n\nComplete success space cadet!! Finished without errors.\n\n')
    
if __name__ == "__main__":
    main()

