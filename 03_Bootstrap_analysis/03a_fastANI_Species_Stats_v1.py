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

import argparse, pathlib
import numpy as np
from collections import defaultdict
from scipy.stats import gaussian_kde
from scipy.signal import find_peaks
from scipy.stats import pearsonr as corr
import matplotlib.pyplot as plt
import seaborn as sns


def parse_fastani(infile, min_ani):

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
    sFRACs = defaultdict(list)
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
            sfrac = vals[1]
            sFRACs[species].append(sfrac)

    return ANIs, sFRACs


def get_peaks(species, anis, sfracs, kde_bandwidth, outpre):

    # run this function per species
    # compute stats and plot peaks kde histogram.
    # list to store data
    # ANI min, q25, median, mean, q50, max
    # Valley_Count, Valley_Widths, Valley_Midpoints.

    # create output directory
    prepath = outpre.split('/')
    if len(prepath) > 1:
        outprx = '/'.join(prepath[:-1])
        outsufx = prepath[-1]
        pathlib.Path(f'{outprx}/peak_plots').mkdir(parents=True, exist_ok=True)
        pathlib.Path(f'{outprx}/nopeak_plots').mkdir(parents=True, exist_ok=True)
        peak_out = f'{outprx}/peak_plots/{outsufx}'
        nopeak_out = f'{outprx}/nopeak_plots/{outsufx}'
    else:
        pathlib.Path('peak_plots').mkdir(parents=True, exist_ok=True)
        pathlib.Path('nopeak_plots').mkdir(parents=True, exist_ok=True)
        peak_out = f'peak_plots/{outpre}'
        nopeak_out = f'nopeak_plots/{outpre}'

    outpeak = f'{peak_out}_{species}.pdf'
    outnopeak = f'{nopeak_out}_{species}.pdf'

    comps = len(anis)
    amin = round(np.min(anis), 4)
    q25 = round(np.quantile(anis, 0.25), 4)
    mean = round(np.mean(anis), 4)
    median = round(np.median(anis), 4)
    q75 = round(np.quantile(anis, 0.75), 4)
    amax = round(np.max(anis), 4)

    kde = gaussian_kde(anis, bw_method=kde_bandwidth)
    x = np.linspace(amin, amax, 1000)
    y = kde(x)
    yn = y * -1 # inverse of y to find the valleys
    yp = find_peaks(y) # returns array (peak indices in x, {properties dict})
    yv = find_peaks(yn) # these are the valley position

    # use the peak indices to grab the values from x
    peaks = sorted(x[yp[0]], reverse=True)
    valleys = sorted(x[yv[0]], reverse=True)

    pc, vc, vws = len(peaks), len(valleys), [] # valley count, valley widths

    if pc > 1: # if more than 1 peak
        # calculate valley widths (ANI distance between peaks)
        for i in range(pc-1):
            p1, p2 = peaks[i], peaks[i+1]
            vws.append(round(p1-p2, 4))
        # append valley midpoint
        vms = [round(i, 4) for i in valleys]

        ovws = ', '.join([str(i) for i in vws])
        ovms = ', '.join([str(i) for i in vms])
        data = [comps, amin, q25, mean, median, q75, amax, vc, ovws, ovms]

        # build the histogram/kde/peak plot

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7,7))

        # ax1 anis, shared frac scatter
        _ = ax1.scatter(anis, sfracs, s=20, alpha=0.2, marker='.')
        ax1.set_xlabel('Average nucleotide identity (%)', fontsize=12)
        ax1.set_ylabel('Shared genome fraction', fontsize=12)

        # ax2 histogram/kde
        bins = np.arange(round(min(anis) - 0.1, 1), round(max(anis), 1) + 0.1, 0.1)
        _ = ax2.hist(anis, bins=bins, density=True)
        _ = ax2.plot(x, y, 'r-', label='Gaussian KDE')
        for mid in vms:
            _ = ax2.axvline(x=mid, color='b', linestyle='--')

        _ = ax2.axvline(x=vms[0], color='b', linestyle='--', label='Midpoint(s)')
        
        ax2.set_xlabel('Average nucleotide identity (%)', fontsize=12)
        ax2.set_ylabel('Density', fontsize=12)
        ax2.legend(loc='best', frameon=False)

        fig.suptitle(' '.join(species.split('_')), fontsize=12)
        fig.set_tight_layout(True)
        plt.savefig(outpeak)
        plt.close()

    else:
        data = [comps, amin, q25, mean, median, q75, amax, vc, 'nan', 'nan']
        vws, vms = None, None

        # build the histogram/kde/ no peak plot

        fig, ax1 = plt.subplots(1, 1, figsize=(7,3.5))

        # ax1 anis, shared frac scatter
        _ = ax1.scatter(anis, sfracs, s=20, alpha=0.2, marker='.')
        ax1.set_xlabel('Average nucleotide identity (%)', fontsize=12)
        ax1.set_ylabel('Shared genome fraction', fontsize=12)

        fig.suptitle(' '.join(species.split('_')), fontsize=12)
        fig.set_tight_layout(True)
        plt.savefig(outnopeak)
        plt.close()

    # data has all stats for the line out to write data table
    # vws is a list of the value widths in descending order
    # vms is a list of the value midpoints in descending order
    return data, vws, vms


def get_metrics(ANIs, sFRACs, min_ani, kde_bd, outpre):

    # compute descriptive stats
    print('\nComputing descriptive stats ...\n')
    # store some values for final plots
    # dictionary for storing the data
    # comps stores number of pairwise genome comparisons
    # valleys stores the number of valleys detected
    # vwidths stores the widths between the detected peaks
    # vmids stores the valley mid points
    # tvm stores the top valley mid points only
    # higap store valley midpoints for species with top valley above 99.8
    # gap stores midpoints for species with top valley in 99-99.8 region
    # nogap stores midpoints for species with top valley below 99.2
    # cc is clonal count when valley is above genomovar range
    # gc is count of species with top valley inside the gap
    # bgc is count of species with top valley below the gap
    data = {
            'comps': [], 'valleys': [], 'vwidths': [], 'vmids': [], 
            'tvm': [], 'gap': [], 'nogap': [], 'higap': [],
            'cc': 0, 'gc': 0, 'bgc': 0
            }


    # setup output file
    data_out = open(f'{outpre}_table.tsv', 'w')
    header = (
            'Species_Name\tComparisons\tMin\tQ25\tMean\tMedian\t'
            'Q50\tMax\tValley_Count\tValley_Widths\tValley_Midpoints\n'
            )
    data_out.write(header)
    removed, higap, gap, logap = [], [], [], []
    for species, anis in ANIs.items():
        # disinclude any species that does not meet the minimum requirements
        if len(anis) < 45:
            removed.append(species)
            continue
        sfracs = sFRACs[species]
        speci, vws, vms = get_peaks(species, anis, sfracs, kde_bd, outpre)
        dline = '\t'.join([str(i) for i in speci])
        data_out.write(f'{species}\t{dline}\n')
        data['comps'].append(speci[0])
        data['valleys'].append(speci[7])
        if vms:
            data['vwidths'].extend(vws)
            data['vmids'].extend(vms)
            top_valley = max(vms)
            data['tvm'].append(top_valley)
            # if top valley above gap range
            if top_valley > 99.8:
                data['higap'].extend(vms)
                data['cc'] += 1
                higap.append(species)
            # if top valley mid point in gap range
            elif top_valley >= 99.2:
                data['gap'].extend(vms)
                data['gc'] += 1
                gap.append(species)
            # else add it to the nogaps
            else:
                data['nogap'].extend(vms)
                data['bgc'] += 1
                logap.append(species)

    data_out.close()

    print(
        f'\n\nThe following {len(removed)} species had fewer than 10 '
        f'pairwise comparisons above the {min_ani} minimum ANI threshold '
        f'and were removed from the analysis:'
        )
    for i in removed:
        print(f'\t\t{i}')

    print(f'\n\n{len(higap)} species with top valley above 99.8% ANI:')
    for i in higap:
        print(f'\t\t{i}')

    print(f'\n\n{len(gap)} species with top valley within [99.2-99.8]% ANI:')
    for i in gap:
        print(f'\t\t{i}')

    print(f'\n\n{len(logap)} species with top valley below 99.2% ANI:')
    for i in logap:
        print(f'\t\t{i}')

    return data


def plot_comps_peaks(comps, peaks, outfile):

    # compute pearson correlation
    pcorr = corr(comps, peaks)
    stats_line = (
        f"Pearson r: {round(pcorr[0], 4)}\n"
        f"p value: {round(pcorr[1], 4)}"
        )

    # Set Colors and markers
    main_color = '#933b41'
    second_color = '#737373'
    vline_color = '#000000'
    color = '#252525'
    marker = '.' #'o'

    # build plot
    g = sns.JointGrid(x=comps, y=peaks)
    print('\nComputing KDEs for marginal plots.')
    # x margin kde plot
    sns.kdeplot(
            x=comps,
            ax=g.ax_marg_x,
            legend=False,
            color=color
            )
    # y margin kde plot
    sns.kdeplot(
            y=peaks,
            ax=g.ax_marg_y,
            legend=False,
            color=color
            )
    # main scatter plot
    g.ax_joint.plot(
        comps,
        peaks,
        marker,
        ms=10,
        alpha=0.25,
        color=color,
        )

    # axis labels and statsline
    g.ax_joint.set_xlabel(
        'Pairwise genome comparisons',
        fontsize=12, y=-0.02
        )
    g.ax_joint.set_ylabel(
        'Peak count',
        fontsize=12, x=-0.02
        )
    g.ax_joint.text(
        0.25, 0.99, stats_line,
        fontsize=10, color=second_color,
        verticalalignment='top', horizontalalignment='right',
        transform=g.ax_joint.transAxes
        )

    # set grid style
    g.ax_joint.yaxis.grid(
        which="major", color='#d9d9d9', linestyle='--', linewidth=1
        )
    g.ax_joint.xaxis.grid(
        which="major", color='#d9d9d9', linestyle='--', linewidth=1
        )
    g.ax_joint.set_axisbelow(True)

    g.ax_joint.set_xscale('log')

    # adjust layout, save, and close
    g.fig.set_figwidth(7)
    g.fig.set_figheight(5)
    g.savefig(outfile)
    plt.close()

    return True


def plot_valley_data(data, outfile):

    # set params
    # comps stores number of pairwise genome comparisons
    # valleys stores the number of valleys detected
    # vwidths stores the widths between the detected peaks
    # vmids stores the valley mid points
    # tvm stores the top valley mid points only
    # higap store valley midpoints for species with top valley above 99.8
    # gap stores midpoints for species with top valley in 99-99.8 region
    # nogap stores midpoints for species with top valley below 99.2
    # cc is clonal count when valley is above genomovar range
    # gc is count of species with top valley inside the gap
    # bgc is count of species with top valley below the gap
    vwidths, vmids = data['vwidths'], data['vmids']
    higap, gap, nogap = data['higap'], data['gap'], data['nogap']
    tvm, cc, gc, bgc = data['tvm'], data['cc'], data['gc'], data['bgc']

    fig, (ax1, ax2, ax3, ax4, ax5, ax6) = plt.subplots(6, 1, figsize=(7,10.5))

    # plot vwidths - valley widths - all species
    bins = np.arange(round(min(vwidths), 1) - 0.1, round(max(vwidths) + 0.1, 1), 0.1)
    _ = ax1.hist(vwidths, bins=bins)
    ax1.set_xlabel('ANI valley widths (%)', fontsize=12)
    ax1.set_ylabel('Count', fontsize=12)
    ax1.text(
        0.99, 0.96, 'All species',
        fontsize=12, color='k',
        verticalalignment='top', horizontalalignment='right',
        transform=ax1.transAxes
        )

    # plot vmids - valley mid points - all species
    bins = np.arange(round(min(vmids), 1) - 0.1 , round(max(vmids) + 0.1, 1), 0.1)
    _ = ax2.hist(vmids, bins=bins)
    #ax2.set_xlabel('ANI valley midpoints (%)', fontsize=12)
    ax2.set_ylabel('Count', fontsize=12)
    ax2.text(
        0.01, 0.96, 'All species',
        fontsize=12, color='k',
        verticalalignment='top', horizontalalignment='left',
        transform=ax2.transAxes
        )

    # plot higap - 1st valley above 99.8 genomovar range
    bins = np.arange(round(min(higap), 1) - 0.1 , round(max(higap) + 0.1, 1), 0.1)
    _ = ax3.hist(higap, bins=bins)
    #ax3.set_xlabel('ANI valley midpoints (%)', fontsize=12)
    ax3.set_ylabel('Count', fontsize=12)
    ax3.text(
        0.01, 0.96, f'{cc} species with top valley above 99.8% ANI',
        fontsize=12, color='k',
        verticalalignment='top', horizontalalignment='left',
        transform=ax3.transAxes
        )

    # plot gap - 1st valley within [99.2-99.8] genomovar range
    bins = np.arange(round(min(gap), 1) - 0.1 , round(max(gap) + 0.1, 1), 0.1)
    _ = ax4.hist(gap, bins=bins)
    ax4.set_ylabel('Count', fontsize=12)
    ax4.text(
        0.01, 0.96, f'{gc} species with top valley within [99.2-99.8]% ANI',
        fontsize=12, color='k',
        verticalalignment='top', horizontalalignment='left',
        transform=ax4.transAxes
        )

    # plot nogap - first valley below 99.2 genomovar range
    bins = np.arange(round(min(nogap), 1) - 0.1 , round(max(nogap) + 0.1, 1), 0.1)
    _ = ax5.hist(nogap, bins=bins)
    ax5.set_ylabel('Count', fontsize=12)
    ax5.text(
        0.01, 0.96, f'{bgc} species with top valley below 99.2% ANI',
        fontsize=12, color='k',
        verticalalignment='top', horizontalalignment='left',
        transform=ax5.transAxes
        )

    # plot tvm
    bins = np.arange(round(min(tvm), 1) - 0.1 , round(max(tvm) + 0.1, 1), 0.1)
    _ = ax6.hist(tvm, bins=bins)
    ax6.set_xlabel('ANI valley midpoints (%)', fontsize=12)
    ax6.set_ylabel('Count', fontsize=12)
    ax6.text(
        0.01, 0.96, f'Top valley midpoints only',
        fontsize=12, color='k',
        verticalalignment='top', horizontalalignment='left',
        transform=ax6.transAxes
        )

    #ax1.legend(loc='best', frameon=False)
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
        help='(OPTIONAL) Specify the minimum ani value (Default = 95).',
        metavar='',
        type=int,
        required=False,
        default=95
        )
    parser.add_argument(
        '-k', '--kde_bandwidth',
        help='(OPTIONAL) Change KDE bandwidth for smoothing (default=0.15).',
        metavar='',
        type=float,
        required=False,
        default=0.15
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\n\nRunning Script...\n')

    # define parameters
    infile= args['input_fastANI_file']
    outpre = args['output_file_prefix']
    min_ani = args['minimum_ani_value']
    kde_bw = args['kde_bandwidth']

    # parse the fastANI file
    print('\nParsing fastANI file ...\n')
    ANIs, sFRACs = parse_fastani(infile, min_ani)

    # compute some metrics, build peak plots per species, write data table
    data = get_metrics(ANIs, sFRACs, min_ani, kde_bw, outpre)

    # plot Comparisons vs Valley_Count
    outfile = f'{outpre}_comp_peak.pdf'
    _ = plot_comps_peaks(data['comps'], data['valleys'], outfile)
    # plot histograms for valley data
    outfile = f'{outpre}_valley_data.pdf'
    _ = plot_valley_data(data, outfile)


    print(f'\n\nComplete success space cadet!! Finished without errors.\n\n')
    
if __name__ == "__main__":
    main()

