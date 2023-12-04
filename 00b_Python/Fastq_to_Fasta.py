#!/usr/bin/env python

''' Converts fastq file to a fasta file.              

Contains read_fastq function whcih can be loaded as a module
into other scripts to return a fastq file as read, seq, qal
objects.

Only works for short read fastq files where each sequence and
quality entry is on a single line.

Does not work with longer read fastq files where sequence
and quality lines may wrap to multiple lines per entry.

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: February 7th, 2020
License :: GNU GPLv3
Copyright 2020 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse

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
    args=vars(parser.parse_args())

    outfile = args['fastq_input_file'].split('.')[0] + '.fasta'

    tracker = {}

    with open(args['fastq_input_file'], 'r') as f, open (outfile, 'w') as o:

        for i, (name, seq, qal) in enumerate(read_fastq(f)):
            if seq >= 10000 and name not in tracker:
                tracker[name] = ''
                o.write(f'>{name}_{i:06}\n{seq}\n')

if __name__ == "__main__":
    main()