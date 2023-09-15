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
    ANI = []
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
            ANI.append(ani)

    return ANIs, ANI


def do_random_sample_test(ANI, ANIs, min_ani, boots, outfile):

    # plot
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(7,10))
    #fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(7,10))

    ###########################################################################
    # plot 1 - distribution of all within species ANI values
    binset = np.arange(94.9, 100.1, 0.1)
    count, bins, ignored = ax1.hist(ANI, binset)
    mn = np.round(np.mean(count), 2)
    amin, amax = np.min(count), np.max(count)
    Umean = 0.5 * (amin + amax)
    ax1.axhline(
            y=np.mean(count), label='Gaussian mean',
            color='#bdbdbd', linestyle='-.', linewidth=2
            )
    ax1.axhline(
            y=Umean, label='Uniform mean',
            color='#636363', linestyle='--', linewidth=2
            )
    ax1.axvline(
            x=min_ani, label='Cutpoint',
            color='k', linestyle='-', linewidth=1
            )
    cutindex = np.where(bins >= min_ani)[0][0]
    cutbins = bins[cutindex:]
    cutcount = count[cutindex:]
    lm = cutbins[np.nonzero(cutcount == np.min(cutcount))[0][-1]] # local minimum
    ax1.axvline(
            x=lm, label='Local minimum',
            color='k', linestyle=':', linewidth=1
            )
    ax1.text(
            0.55, 0.96, f'ANI distribution (all data)',
            fontsize=12, color='k',
            verticalalignment='top', horizontalalignment='center',
            transform=ax1.transAxes
            )
    ax1.legend(loc='best')
    ax1.set_ylabel('Count', fontsize=12)
    ###########################################################################

    binset = np.arange(min_ani, 100.1, 0.1)
    boot_sample, kdes = bootstraps(ANIs, min_ani, boots)

    ###########################################################################
    # plot 2 10% random sample
    count, bins, ignored = ax2.hist(boot_sample, binset, density=True)

    kde = gaussian_kde(boot_sample, bw_method=0.15)
    x = np.linspace(min_ani, 100, 1000)
    y = kde(x)
    yn = y * -1 # inverse of y to find the valleys
    yp = find_peaks(y) # returns array (peak indices in x, {properties dict})
    yv = find_peaks(yn) # these are the valley position
    # use the peak indices to grab the values from x
    peaks = sorted(x[yp[0]], reverse=True)
    valleys = sorted(x[yv[0]], reverse=True)

    _ = ax2.plot(x, y, 'r-', label='Gaussian KDE')

    ax2.axvline(
            x=peaks[0], label='Local maximums',
            color='#bdbdbd', linestyle='--', linewidth=1
            )
    ax2.axvline(
            x=valleys[0], label='Local minimums',
            color='#636363', linestyle='-.', linewidth=1
            )
    if len(peaks) > 1:
        for pk in peaks[1:]:
            ax2.axvline(
                x=pk, color='#bdbdbd', linestyle='--', linewidth=1
                )
    if len(valleys) > 1:
        for val in valleys[1:]:
            ax2.axvline(
                x=val, color='#636363', linestyle='-.', linewidth=1
                )
    ax2.text(
        0.6, 0.96, f'ANI distribution (bootstrap sample)',
        fontsize=12, color='k',
        verticalalignment='top', horizontalalignment='center',
        transform=ax2.transAxes
        )
    ax2.legend(loc='best')
    ax2.set_ylabel('Density', fontsize=12)
    ###########################################################################

    ###########################################################################
    # plot 3 distribution of bootstrapped confidence interval
    dfy = pd.DataFrame(kdes['y'])
    x, y1, y2, y3 = kdes['x'][0], dfy.mean(), dfy.quantile(0.025), dfy.quantile(0.975)
    _ = ax3.plot(x, y1, 'b--', alpha=1)
    _ = ax3.plot(x, y2, 'b:', alpha=0.5)
    _ = ax3.plot(x, y3, 'b:', alpha=0.5)
    ax3.fill_between(x, y2, y3, alpha=0.2)

    # itertools.chain.from_iterable() flattens a list of lists to single list
    all_peaks = list(itertools.chain.from_iterable(kdes['peaks']))
    all_valleys = list(itertools.chain.from_iterable(kdes['valleys']))
    for pk in all_peaks:
            ax3.axvline(
                x=pk, color='#bdbdbd', linestyle='--', linewidth=1, alpha=0.1
                )
    for val in all_valleys:
            ax3.axvline(
                x=val, color='#636363', linestyle='-.', linewidth=1, alpha=0.1
                )
    ax3.text(
        0.01, 0.96, f'Bootstrapped\nconfidence\nintervals',
        fontsize=12, color='k',
        verticalalignment='top', horizontalalignment='left',
        transform=ax3.transAxes
        )
    ax3.set_ylabel('Density', fontsize=12)
    ###########################################################################

    ###########################################################################
    # plot 4 distribution of bootstrapped local minimums and maximums
    _ = ax4.hist(all_peaks, binset, color='#bdbdbd', label='Local maximums')
    _ = ax4.hist(all_valleys, binset, color='#636363', label='Local minimums')

    ax4.set_ylabel('Count', fontsize=12)
    ax4.set_xlabel('Average nucleotide identity (%)', fontsize=12)
    ax4.legend(loc='best')
    ###########################################################################

    fig.set_tight_layout(True)
    plt.savefig(outfile)
    plt.close()

    return True


def bootstraps(ANIs, min_ani, boots):

    # Random sampling
    np.random.seed(int.from_bytes(os.urandom(4), byteorder='little'))

    # step 1: Select only species with 100 pairwise comparisons above 96.5
    select = [] # store selected species
    for species, ani in ANIs.items():
        bani = [i for i in ani if i >= min_ani]
        if len(bani) >= 100 and np.mean(bani) < 99.5:
            select.append(bani)
    print(f'\nNon-clonal species with ≥ 100 pairwise comparisons: {len(select)}')

    # step 2: Random sample with replacement 100 points from each selected species
    kdes = {'x': [], 'y': [], 'peaks': [], 'valleys': []}
    timer = 0
    for i in range(boots):
        boot_sample = []
        for species in select:
            s = np.random.choice(species, 100, replace=True)
            boot_sample.extend(s)
        kde = gaussian_kde(boot_sample, bw_method=0.15)
        x = np.linspace(min_ani, 100, 1000)
        y = kde(x)
        yn = y * -1 # inverse of y to find the valleys
        yp = find_peaks(y) # returns array (peak indices in x, {properties dict})
        yv = find_peaks(yn) # these are the valley position
        # use the peak indices to grab the values from x
        kdes['x'].append(x)
        kdes['y'].append(y)
        kdes['peaks'].append(sorted(x[yp[0]], reverse=True))
        kdes['valleys'].append(sorted(x[yv[0]], reverse=True))

        timer += 1
        if timer % 20 == 0:
            print(f'\t\t... Bootstrap {timer}')

    return boot_sample, kdes


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
    parser.add_argument(
        '-boots', '--bootstrap_count',
        help='(OPTIONAL) Number of bootstrap iterations (Default = 1000).',
        metavar='',
        type=int,
        required=False,
        default=1000
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\n\nRunning Script...\n')

    # define parameters
    infile= args['input_fastANI_file']
    outpre = args['output_file_prefix']
    min_ani = args['minimum_ani_value']
    boots = args['bootstrap_count']

    # parse the fastANI file
    print('\nParsing fastANI file ...\n')
    ANIs, ANI = parse_fastani(infile)

    # plot histograms for valley data
    print('\nRunning bootstrapped tests ...\n')
    outfile = f'{outpre}_Bootstrap_test.pdf'
    _ = do_random_sample_test(ANI, ANIs, min_ani, boots, outfile)

    print(f'\n\nComplete success space cadet!! Finished without errors.\n\n')
    
if __name__ == "__main__":
    main()

