 #!/usr/bin/env python

'''Reads all vs all fastANI output and computes descriptive statistics

fastANI output includes the full path to query and reference genomes.

This script expects the following file path naming scheme:

    genomes/Species_Name/genome.fasta

Where the script cuts out "Species_Name" to group the genomes.

The only input is the all vs all fastANI file and the output is a tsv
file containing the following columns:

Species_Name, Comparisons, ANI min, q25, median, mean, q50, max,
Valley_Count, Valley_Widths, Valley_Midpoints.

Comparisons is the number of genome pair combinations per species.
ANI min, q25, median, mean, q50, max are self explanatory.
Valley_Count is the number of valleys detected in the Gaussian KDE fit to
a histogram of the ANI distribution for each species.
Valley_Widths is the distance in ANI between the peaks.
Valley_Midpoint is the ANI value at the midpoint between peaks.
Valleys_Widths and Valley_Midpoints are only reported if more than 1
peak is detected.

Creates directory: peak_plots
Writes out histogram/kde/peak plots for each species.

Writes Valley_Count vs genome count plot.

Writes valley midpoint histogram.

bins size 0.1 ANI.

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: Sept. 2023
License :: GNU GPLv3
Copyright 2023 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse, pathlib, os, itertools
from collections import defaultdict
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import seaborn as sns


def parse_fastani(infile):

    min_ani = 95 # only interesting in above the species level.

    # Parses fastANI output file into a dict of {species: [ANI scores]}
    # Removes self matched scores and retains only one ANI score per
    # genome pair by keeping the ANI score with the higher total_frags
    # because each genome can be used as query and reference.

    # first dict to store all ANI values and filter for higher total_frags
    # dict of {'-'.join(sorted([qry,ref])): total_frags}
    data = {}

    with open(infile, 'r') as file:
        for line in file:
            X = line.rstrip().split('\t')
            # get genome file names and remove self matches
            qry = X[0]
            ref = X[1]
            if qry == ref: continue # remove self matches
            gpair = ';;;'.join(sorted([qry, ref]))
            # get ANI
            ani = float(X[2])
            # shared fraction used to retain only 1 ANI value per genome
            # pair using the ANI value from the largest genome as the
            # reference genome. Largest as defined by total_frags.
            shared_frags = int(X[3])
            total_frags = int(X[4])
            shared_frac = shared_frags / total_frags

            if ani >= min_ani:

                if gpair in data:
                    old_frags = data[gpair][2]
                    if total_frags > old_frags:
                        data[gpair] = [ani, shared_frac, total_frags]
                else:
                    data[gpair] = [ani, shared_frac, total_frags]

    # iterate through filter genome/ANI data and store by species
    # create dict to store {species: [ANI scores]}
    ANIs = defaultdict(list)
    for gpair, vals in data.items():
            # get species names
            X = gpair.split(';;;')
            species1 = X[0].split('/')[1]
            species2 = X[1].split('/')[1]
            if species1 == species2:
                species = species1
            else:
                species = '-'.join(sorted([species1, species2]))
            ani = vals[0]
            ANIs[species].append(ani)

    return ANIs


def plot_species_bin_density(ANIs, min_ani, outpre):

    binset = np.round(np.arange(min_ani-0.1, 100.1, 0.1), 1)

    # hiANI is species mean ANI ≥ 99.5
    # loANI_A is species has a local minimum in range [99.2-99.8]
    # loANI_B is species has a local maximum in range [99.2-99.8]
    # loANI_C is species without local min or max in [99.2-99.8]
    hiANI, loANI_A, loANI_B,loANI_C = {}, {}, {}, {}

    for species, anis in ANIs.items():
        # skip species if fewer than 45 pairwise comparisons above 95% ANI.
        # this comes from the dict created from the parse_fastani function.
        if len(anis) < 45: continue
        # filter ani set to min_ani
        fanis = [i for i in anis if i >= min_ani]
        # create histogram and get the counts
        count, edges = np.histogram(fanis, binset)
        # sort by local maximum in "gap" range
        kde = gaussian_kde(fanis, bw_method=0.15)
        x = np.linspace(min_ani, 100, 1000)
        y = kde(x)
        yn = y * -1 # inverse of y to find the valleys
        yp = find_peaks(y) # returns array (peak indices in x, {properties dict})
        yv = find_peaks(yn) # these are the valley position
        # use the peak indices to grab the values from x
        peaks = sorted(x[yp[0]], reverse=True)
        valleys = sorted(x[yv[0]], reverse=True)
        # select any valleys in the "gap" range.
        vly = [i for i in valleys if i <= 99.8 and i >= 99.2]
        pks = [i for i in peaks if i <= 99.8 and i >= 99.2]
        # sort by mean ANI value
        mean_ani = np.mean(fanis)
        # replace species name '_' with space
        name = species.replace('_', ' ')
        # Classify hiANI or "clonal" species
        if mean_ani >= 99.5:
            # normalize counts to frequency and convert to percent
            hiANI[name] = (count/len(fanis)) * 100
        # Classify loANI or heterogeneous species
        # Classify loANI into:
        # A: has peak in gap region
        # or B: no peak in gap region
        else:
            # loANI A will have a valley in vly
            if len(vly) >= 1:
                # normalize counts to frequency and convert to percent
                loANI_A[name] = (count/len(fanis)) * 100
            else: # loANI B will not have any valleys in vly but will have peak
                if len(pks) >= 1:
                    # normalize counts to frequency and convert to percent
                    loANI_B[name] = (count/len(fanis)) * 100
                else:
                    # no peak or valley or in range [99.2-99.8]
                    loANI_C[name] = (count/len(fanis)) * 100
        
    # hiANI group plot
    # create df
    df = pd.DataFrame(hiANI, index=binset[1:])
    print(f'\t\tClonal species: {len(df.columns)}')
    g = sns.clustermap(
                        df.T, figsize=(36,48),
                        xticklabels=True,
                        yticklabels=True,
                        col_cluster=False
                        )

    # print out clustered figure labels
    #texts = [t.get_text()  for t in g.ax_heatmap.get_yticklabels()]
    #for i in texts: print(i)
    outfile = f'{outpre}_hiANI_histclust.pdf'
    g.savefig(outfile)
    plt.close()

    # loANI group A plot
    # create df
    df = pd.DataFrame(loANI_A, index=binset[1:])
    print(f'\t\tNon-clonal pattern A species: {len(df.columns)}')
    g = sns.clustermap(
                        df.T, figsize=(36,48),
                        xticklabels=True,
                        yticklabels=True,
                        col_cluster=False
                        )
    outfile = f'{outpre}_loANI_A_histclust.pdf'
    g.savefig(outfile)
    plt.close()

    # loANI group B plot
    # create df
    df = pd.DataFrame(loANI_B, index=binset[1:])
    print(f'\t\tNon-clonal pattern B species: {len(df.columns)}')
    g = sns.clustermap(
                        df.T, figsize=(36,48),
                        xticklabels=True,
                        yticklabels=True,
                        col_cluster=False
                        )
    outfile = f'{outpre}_loANI_B_histclust.pdf'
    g.savefig(outfile)
    plt.close()

    # loANI group C plot
    # create df
    df = pd.DataFrame(loANI_C, index=binset[1:])
    print(f'\t\tNon-clonal pattern C species: {len(df.columns)}')
    g = sns.clustermap(
                        df.T, figsize=(36,48),
                        xticklabels=True,
                        yticklabels=True,
                        col_cluster=False
                        )
    outfile = f'{outpre}_loANI_C_histclust.pdf'
    g.savefig(outfile)
    plt.close()

    return True


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--input_fastANI_file',
        help='Please specify the input all vs all fastANI file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--output_file_prefix',
        help='Please specify prefix for output file names!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-min_ani', '--minimum_ani_value',
        help='(OPTIONAL) Minimum ani value for bootstraps (Default = 96.5).',
        metavar='',
        type=float,
        required=False,
        default=96.5
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\n\nRunning Script...\n')

    # define parameters
    infile= args['input_fastANI_file']
    outpre = args['output_file_prefix']
    min_ani = args['minimum_ani_value']

    # parse the fastANI file
    print('\nParsing fastANI file ...\n')
    ANIs = parse_fastani(infile)

    # plot clustermap of bin densities
    print('\nBuilding bin density clustermaps ...\n')
    _ = plot_species_bin_density(ANIs, min_ani, outpre)

    print(f'\n\nComplete success space cadet!! Finished without errors.\n\n')
    
if __name__ == "__main__":
    main()

