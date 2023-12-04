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

import argparse, shutil
from pathlib import Path

def gather_data(infile, outpre, directory):
    # Creates all vs all file for clinical vs environmental classifications

    # initialize lists to store species names
    clinical = []
    environmental = []

    # read through the tsv file and populate lists
    with open(infile, 'r') as file:
        header = file.readline()
        for line in file:
            X = line.rstrip().split('\t')
            name = X[0].split('.')[0].replace(" ", "_")
            clss1 = X[1]
            clss2 = X[2] if len(X) > 2 else 'na'

            if clss1 == 'environmental':
                environmental.append(name)
                print(f'{name} assigned to environmental.')
            elif clss1 == 'clinical' and clss2 == 'human':
                clinical.append(name)
                print(f'{name} assigned to clinical.')
            else:
                print(f'{name} not assigned.')
    
    # concatenate clinical ani results
    with open(f'{outpre}_clinical.ani', 'wb') as outfile:
        for name in clinical:
            ani_file = f'{directory}/{name}.ani'
            if Path(ani_file).is_file():
                with open(ani_file, 'rb') as fb:
                    shutil.copyfileobj(fb, outfile)

    # concatenate environmental ani results
    with open(f'{outpre}_environmental.ani', 'wb') as outfile:
        for name in environmental:
            ani_file = f'{directory}/{name}.ani'
            if Path(ani_file).is_file():
                with open(ani_file, 'rb') as fb:
                    shutil.copyfileobj(fb, outfile)


    return True


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--input_file',
        help='Please specify the input classification tsv file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--output_prefix',
        help='Please specify the output file prefix!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-d', '--ANI_results_directory',
        help='Please specify the output file prefix!',
        metavar='',
        type=str,
        required=True
        )

    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\n\nRunning Script...\n')

    # define parameters
    infile = args['input_file']
    outpre = args['output_prefix']
    directory = args['ANI_results_directory']


    # read in the data
    _ = gather_data(infile, outpre, directory)

    print(f'\n\nComplete success space cadet!! Finished without errors.\n\n')

if __name__ == "__main__":
    main()
