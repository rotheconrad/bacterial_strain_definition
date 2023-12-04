#!/usr/bin/env python

''' Converts reads in a fastq file to a separate fasta file for each read.              

Contains read_fastq function which can be loaded as a module
into other scripts to return a fastq file as read, seq, qal
objects.

Only works for fastq files where each sequence and quality entry is on a
single line.

Will not work if sequences are on more than one line.

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: September, 2022
License :: GNU GPLv3
Copyright 2022 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse
from pathlib import Path

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

def main():
        # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--fastq_input_file',
        help='Please specify the fastq file to read!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--fasta_output_file',
        help='Please specify the fasta output file name!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-d', '--output_directory',
        help='Please specify the output directory for single fasta!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-m', '--minimum_read_length',
        help='(OPTIONAL) Specify the minimum read length (Default = 10000)!',
        metavar='',
        type=int,
        required=False,
        default=10000
        )
    args=vars(parser.parse_args())

    infile = args['fastq_input_file']
    outfile = args['fasta_output_file']
    outdir = args['output_directory']
    lmin = args['minimum_read_length']

    # make outdir
    Path(outdir).mkdir(parents=True, exist_ok=True)
    # keep track of read names and remove duplicates
    tracker = {}
    count = 1
    with open(infile, 'r') as f, open(outfile, 'w') as o:
        for name, seq, qal in read_fastq(f):
            if len(seq) >= lmin:
                if name in tracker:
                    old_seq = tracker[name]
                    if seq == old_seq:
                        continue
                    else:
                        name = f'{name}-{count}'

                tracker[name] = seq
                o.write(f'>{name}\n{seq}\n')
                singleout = f'{outdir}/read_{count:06}.fasta'
                with open(singleout, 'w') as so:
                    so.write(f'>{name}\n{seq}\n')
                
            count += 1

if __name__ == "__main__":






    main()