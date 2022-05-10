#!/usr/bin/env python

'''Parse NCBI Assembly Summary File.

This scripts takes an NCBI assembly summary file such as:

ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt
or
ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt

This script outputs files for counts of each species at each assembly
level and files to download genomes for all species at each assembly
level. Filters for prokaryotes.

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: Apr 2022
License :: GNU GPLv3
Copyright 2022 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse
from collections import defaultdict


def parse_ncbi_assembly_summary(infile, outpre, n):

    # create dictionary to store the data
    data = {
            "Complete": defaultdict(list),
            "Chromosome": defaultdict(list),
            "Scaffold": defaultdict(list),
            "Contig": defaultdict(list),
            }

    # read file, populate data dict
    with open(infile, 'r') as file:
        for line in file:
            # skip header
            if line.startswith('#'): continue
            # split by tabs
            X = line.rstrip().split('\t')
            # select species name but not strain
            organism = ' '.join(X[7].split(' ')[:2])
            # get the version number to filter for latest
            version = X[10]
            # select assembly level
            level = X[11].split(' ')[0]
            # get ftp address for genome fasta file
            ftp = X[19] + '/' + X[19].split('/')[-1]

            # select the latest version
            if version == 'latest':
                data[level][organism].append(ftp)

    # Read through data dict count genomes at each level for each species
    # and write ftp file for each level:
    for levels, organisms in data.items():
        f1 = f'{outpre}_{levels}_counts.tsv'
        f2 = f'{outpre}_{levels}_ftps.sh'
        with open(f1, 'w') as f1out, open(f2, 'w') as f2out:
            for species, ftps in organisms.items():
                count = len(ftps)
                name = species.replace(" ", "_").replace(".", "")
                print(f'{levels}\t{species}\t{count}')
                f1out.write(f'{species}\t{count}\n')
                if count >= n:
                    for i in ftps:
                        f2out.write(f'{name}\t{i}_genomic.fna.gz\n')

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
        '-p', '--output_prefix',
        help='Please specify the output file prefix!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-n', '--genome_count_filter',
        help='Optional: Genome count filter (default=100)!',
        metavar='',
        type=int,
        required=False,
        default=100
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\n\nRunning Script...\n')

    infile = args['input_file']
    outpre = args['output_prefix']
    n = args['genome_count_filter']

    _ = parse_ncbi_assembly_summary(infile, outpre, n)

    print(f'\n\nComplete success space cadet!! Finished without errors.\n\n')

if __name__ == "__main__":
    main()
