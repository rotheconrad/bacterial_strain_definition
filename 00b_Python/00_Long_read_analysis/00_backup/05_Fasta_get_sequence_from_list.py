#!/usr/bin/env python

''' Retreive sequences from fasta file from input list of seq names

Input is a list of seq names, 1 name per line include the '>' and a
fasta file.

Returns a new fasta file with only the sequences from the list.

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

import argparse, subprocess

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

def parse_input_list(input_list):

    seq_names = {}

    with open(input_list, 'r') as file:
        for line in file:
            seq_names[line.rstrip()] = ''

    return seq_names

def return_matching_fasta(input_fasta, seq_names, output_file):

    with open(input_fasta, 'r') as file, open(output_file, 'w') as outfile:
        for name, seq in read_fasta(file):
            if name in seq_names:
                outfile.write(f'{name}\n{seq}\n')

    return True

def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-l', '--input_list',
        help='Please specify the file path to the list of seq names!',
        metavar=':',
        type=str,
        required=True
        )
    parser.add_argument(
        '-f', '--input_fasta',
        help='Please specify the fasta file path!',
        metavar=':',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--output_file',
        help='Please specify the output file path!',
        metavar=':',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\n\nRunning Script...\n')

    # define parameters
    input_list = args['input_list']
    input_fasta = args['input_fasta']
    output_file = args['output_file']

    # Get the list of seq names
    seq_names = parse_input_list(input_list)

    # Find the matching fasta entries and write to file
    _ = return_matching_fasta(input_fasta, seq_names, output_file)

    print(f'\n\nComplete success space cadet!! Finished without errors.\n\n')
    
if __name__ == "__main__":
    main()

