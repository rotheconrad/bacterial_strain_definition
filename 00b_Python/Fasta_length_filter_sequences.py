#!/usr/bin/env python

''' Filters fasta sequences by length.

Writes new fasta file with all sequences >= minimum length.

The length is in base pairs/characters.

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: June 21st, 2019
License :: GNU GPLv3
Copyright 2019 Roth Conrad
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

def Fasta_rename_sequences(infile, outfile, minlength):

    passing_sequences = 0
    filtered_sequences = 0

    with open(infile, 'r+') as f, open(outfile, 'w') as o:
        for name, seq in read_fasta(f):
            if len(seq) >= minlength:
                passing_sequences += 1
                o.write(f'{name}\n{seq}\n')
            else:
                filtered_sequences += 1

    print(
        f'\n\t\tSequences written to new file: {passing_sequences}\n'
        f'\t\tSequences below length filter: {filtered_sequences}'
        )

    return True

def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--input_file',
        help='Please specify the input fasta file!',
        metavar=':',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--output_file',
        help='Please specify the output file name!',
        metavar=':',
        type=str,
        required=True
        )
    parser.add_argument(
        '-m', '--minimun_sequence_length',
        help='(OPTIONAL) Specify minimum sequence length (default = 10000)!',
        metavar=':',
        type=str,
        required=False,
        default=10000
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\n\nReading fasta file and filtering sequences by length ...\n')

    # Run the renaming function:
    Fasta_length_filter_sequences(
                    args['input_file'],
                    args['output_file'],
                    args['minimun_sequence_length']
                    )

    print(f'\n\nComplete success space cadet!! Finished without errors.\n\n')
    
if __name__ == "__main__":
    main()

