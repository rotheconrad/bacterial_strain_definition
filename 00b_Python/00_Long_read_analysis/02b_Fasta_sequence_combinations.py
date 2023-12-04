#!/usr/bin/env python

''' Returns n choose 2 number of combinations
where n is the number of fasta sequences in
a fasta file.

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: Sept. 2022
License :: GNU GPLv3
Copyright 2022 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse
from math import factorial


def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))


def blocks(files, size=65536):
    while True:
        b = files.read(size)
        if not b: break
        yield b


def count_fasta_entries(infile):
    with open(infile, 'r') as f:
        n = sum(bl.count('>') for bl in blocks(f))

    return n


def compute_combinations(infile):

    # Read fasta file and return sequence count
    n = count_fasta_entries(infile)
    r = 2
    # use formula n choose r = n!/(r!(n-r)!)
    # calculate number of pairwise combinations
    c = factorial(n) / (factorial(r) * factorial(n-r))

    name = infile.split('/')[-1].split('.')[0]
    print(f'{name} Reads: {n} Combinations: {c}')

    return True


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--input_file',
        help='Please specify the paired read interleaved fasta file!',
        metavar='',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())


    infile = args['input_file']
    _ = compute_combinations(infile)

if __name__ == "__main__":
    main()
