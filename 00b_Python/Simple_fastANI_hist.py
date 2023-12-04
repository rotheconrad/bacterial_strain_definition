 #!/usr/bin/env python

'''Returns a histogram of percent sequence similarity from a PAF file.

Intended for reads mapped to a reference MAG or genome.
Generates two histograms of the percent sequence similarity distribution
of the mapped reads with x-axis 70-100 and 98-100. Generates histogram of
read length distribution.

Sequence similarity based on minimap2 manual:
https://lh3.github.io/minimap2/minimap2.html
https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity

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
import seaborn as sns
sns.set_theme()

def parse_fastani(ani_file, xmin, xmax):
    """Reads in the tsv ALLvsAll ANI """

    print("\nReading data.")
    name_dict = {}
    anis = []

    with open(ani_file, 'r') as f:
        for l in f:
            # split each line by tabs
            X = l.rstrip().split('\t')
            # get the genome pair names
            qry_genome = X[0]
            ref_genome = X[1]
            # total frac column. Keep larger genome as reference
            tfrac = int(X[4]) 
            # remove self matchs
            if qry_genome == ref_genome: continue
            # compute metrics
            ani = float(X[2])
            # get the shared genome fraction
            shared = float(X[3])
            total = float(X[4])
            ratio = shared / total
            # sort genome file names and combine
            names = [qry_genome, ref_genome]
            names.sort()
            gname = '-'.join(names)
            # record values from only one genome pair.
            # Keep larger genome as reference.
            if tfrac > name_dict.get(gname[0], 0):
                # add genome pair to dict
                name_dict[gname] = [tfrac, ani, ratio]

    # write data to arrays
    for gpair, metrics in name_dict.items():
        ani = metrics[1]
        if ani > xmin:
            anis.append(ani)

    return anis

def build_pid_plots(anis, xmin, xmax, bins, title, outfile):

    fntsz = 12 # font size
    color = '#639efc' #'#bdbdbd' # bar color

    fig, ax = plt.subplots(figsize=(7,5))

    _ = ax.hist(anis, bins=bins, color=color, alpha=1)
    ax.set_ylabel('Count of comparisons', fontsize=12)
    ax.set_xlabel('FastANI estimate')
    ax.set_title(title, fontsize=fntsz, y=1.1)
    ax.set_xlim(xmin-0.2, xmax+0.2)
    ax.xaxis.grid(visible=False)

    fig.set_tight_layout(True)
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
        '-i', '--input_file',
        help='Please specify the fastANI file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--output_file',
        help='Please specify the output file name (use .pdf)!',
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
    parser.add_argument(
        '-b', '--number_of_bins',
        help='(OPTIONAL) Specify the number of bins (Default = 25).',
        metavar='',
        type=str,
        required=False,
        default=30
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\n\nRunning Script...\n')

    # define parameters
    infile= args['input_file']
    outfile = args['output_file']
    xmin = args['xaxis_minimum']
    xmax = args['xaxis_maximum']
    bins = args['number_of_bins']


    # parse the sam file and compute percent sequence similarity (pID)
    print('\nParsing fastANI file ...\n')
    anis = parse_fastani(infile, xmin, xmax)

    print('\nPlotting histogram ...\n')
    n = len(anis)
    title = f'{n} within species pairwise comparisons\nfrom 330 species closed genomes'
    _ = build_pid_plots(anis, xmin, xmax, bins, title, outfile)

    print(f'\n\nComplete success space cadet!! Finished without errors.\n\n')
    
if __name__ == "__main__":
    main()

