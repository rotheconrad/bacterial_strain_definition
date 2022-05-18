 #!/usr/bin/env python

'''reports fraction of ANI values between a range above a cutoff.



-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: April 2022
License :: GNU GPLv3
Copyright 2022 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse


def gather_data(ani, xmin, xmax):
    """Reads in the tsv ALLvsAll ANI """

    print("\nReading data.")
    name_dict = {}
    xs = []

    with open(ani, 'r') as f:
        for l in f:
            X = l.rstrip().split('\t')

            # get the genome pair names
            qry_genome = X[0]
            ref_genome = X[1]
            # remove self matchs
            if qry_genome == ref_genome: continue
            # sort genome file names and combine
            names = [qry_genome, ref_genome]
            names.sort()
            gname = '-'.join(names)
            # record values from only one genome pair
            if gname in name_dict: continue
            # add genome pair to dict
            name_dict[gname] = ''

            ani = float(X[2])

            xs.append(ani)

    inrange = len([ i for i in xs if i >= xmin and i <= xmax ])
    total = len(xs)
    total95 = len([i for i in xs if i >= 95])
    total96 = len([i for i in xs if i >= 96])
    total97 = len([i for i in xs if i >= 97])
    total98 = len([i for i in xs if i >= 98])
    total99 = len([i for i in xs if i >= 99])

    print(
        '\nGenome pair counts:\n'
        f'\t(A) Genome pairs in range [{xmin}, {xmax}]: {inrange}\n'
        f'\t(B) Total genome pairs: {total}\n'
        f'\t(C) Genome pairs >= 95% ANI: {total95}\n'
        f'\t(D) Genome pairs >= 96% ANI: {total96}\n'
        f'\t(E) Genome pairs >= 97% ANI: {total97}\n'
        f'\t(F) Genome pairs >= 98% ANI: {total98}\n'
        f'\t(G) Genome pairs >= 99% ANI: {total99}\n'
        '\nFraction of A in B-G:\n'
        f'\t(A) / (B) = {inrange/total:.4f}\n'
        f'\t(A) / (C) = {inrange/total95:.4f}\n'
        f'\t(A) / (D) = {inrange/total96:.4f}\n'
        f'\t(A) / (E) = {inrange/total97:.4f}\n'
        f'\t(A) / (F) = {inrange/total98:.4f}\n'
        f'\t(A) / (G) = {inrange/total99:.4f}\n'
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
        help='Please specify the input file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-xmin', '--xaxis_minimum',
        help='OPTIONAL: Minimum value to plot on x-axis. (Default=95.0)',
        metavar='',
        type=float,
        default=95.0,
        required=False
        )
    parser.add_argument(
        '-xmax', '--xaxis_maximum',
        help='OPTIONAL: Maximum value to plot on x-axis. (Default=100.0)',
        metavar='',
        type=float,
        default=100.0,
        required=False
        )

    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\n\nRunning Script...\n')

    # define parameters
    infile = args['input_file']
    xmin = args['xaxis_minimum']
    xmax = args['xaxis_maximum']


    # read in the data
    _ = gather_data(infile, xmin, xmax)

    print(f'\n\nComplete success space cadet!! Finished without errors.\n\n')

if __name__ == "__main__":
    main()
