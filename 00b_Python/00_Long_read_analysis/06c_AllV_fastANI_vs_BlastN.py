 #!/usr/bin/env python

''' Compare all vs all results from fastANI to BlastN

Reads and filters equivalent all vs all fastANI and BlastN tabular 
output files. Calculates the differences between the scores of pairwise
comparisons.

Requires a read pair to file to key. FastANI was built for genomes, so
each sequence to compare needs its own fasta file and the output file
contains file names for query and reference whereas BlastN accepts and 
compares multiple sequences in a fasta file and the default read names
are used for the query and reference. We need them to match. I created
a two column text file separated with a single white space. The file name
(and the path if it is in the fastANI output file) in column 1 and the
corresponding read name (this should be the fasta defline) in column 2.
This script expects column two to have the '>' character of fasta files.

ex:
directory/path/read_014437.fasta >m64078_201009_211857

Outputs tsv file with ANI, pID and difference for each
pairwise comparison found in both results.

Outputs results from fastANI only.

Outputs results from BlastN only.

Output header key:
    - Read pair - name of the matched read pair
    - x diff: ANI - pID
    - y diff: shared fraction - aligned fraction
    - ANI: Score from fastANI
    - pID: Score from BlastN
    - sfrac: shared fraction
    - afrac: aligned fraction

    * pID is the percent sequence alignment from BlastN
    * shared fraction = shared fragments / total fragments
    * aligned fraction = alignment length / query length

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

import argparse, math
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns; sns.set(style="white", color_codes=True)

class TickRedrawer(matplotlib.artist.Artist):
    #https://stackoverflow.com/questions/19677963/
    #matplotlib-keep-grid-lines-behind-the-graph-but-the-y-and-x-axis-above
    """Artist to redraw ticks."""

    __name__ = "ticks"

    zorder = 10

    @matplotlib.artist.allow_rasterization
    def draw(self, renderer: matplotlib.backend_bases.RendererBase) -> None:
        """Draw the ticks."""
        if not self.get_visible():
            self.stale = False
            return

        renderer.open_group(self.__name__, gid=self.get_gid())

        for axis in (self.axes.xaxis, self.axes.yaxis):
            loc_min, loc_max = axis.get_view_interval()

            for tick in axis.get_major_ticks() + axis.get_minor_ticks():
                if tick.get_visible() and loc_min <= tick.get_loc() <= loc_max:
                    for artist in (tick.tick1line, tick.tick2line):
                        artist.draw(renderer)

        renderer.close_group(self.__name__)
        self.stale = False


def file_to_read_key(key_file):
    # creates a dict from the two column input file
    # {file-name: read-name}
    # returns dict

    fast_key = {}

    with open(key_file, 'r') as file:
        for line in file:
            X = line.rstrip().split(' ')
            fname = X[0]
            rname = X[1][1:]
            fast_key[fname] = rname

    return fast_key


def parse_fastani(infile, min_ani, min_frac, fast_key):

    # Reads fastani file, filters results, returns dict

    # {read-pair: [ani, shared_frac]}
    fast_dict = {}

    with open(infile, 'r') as file:
        for line in file:
            X = line.rstrip().split('\t')
            qry = fast_key[X[0]]
            ref = fast_key[X[1]]
            # check for self matches and discard
            if qry == ref: continue
            ani = float(X[2])
            shared_frags = int(X[3])
            total_frags = int(X[4])
            shared_frac = shared_frags / total_frags

            if shared_frac >= min_frac and ani >= min_ani:
                # check for reciprocal matches and discard lower score
                sorted_pair = '-'.join(sorted([qry, ref]))
                if sorted_pair in fast_dict:
                    old_ani = fast_dict[sorted_pair][0]
                    if ani > old_ani:
                        fast_dict[sorted_pair] = [ani, shared_frac]
                else:
                    fast_dict[sorted_pair] = [ani, shared_frac]

    return fast_dict


def parse_blastn(infile, min_pid, min_frac):

    # Reads tabular blastn file, filters results, returns dict

    # {read-pair: [pid, aligned_frac]}
    blast_dict = {}

    with open(infile, 'r') as file:
        for line in file:
            X = line.rstrip().split('\t')
            qname = X[0]
            rname = X[1]
            # check for self matches and discard
            if qname == rname: continue
            pid = float(X[2])
            bitscore = float(X[11])
            alen = int(X[3]) # alignment length
            qlen = int(X[12]) # reference sequence length
            aligned_frac = alen / qlen

            if aligned_frac >= min_frac and pid >= min_pid:
                # check for reciprocal matches and discard lower score
                sorted_pair = '-'.join(sorted([qname, rname]))
                if sorted_pair in blast_dict:
                    old_bitscore = blast_dict[sorted_pair][2]
                    if bitscore > old_bitscore:
                        blast_dict[sorted_pair] = [pid, aligned_frac, bitscore]
                else:
                    blast_dict[sorted_pair] = [pid, aligned_frac, bitscore]

    return blast_dict


def compute_difference(fastANI, BlastN, outpre):

    '''
    Read through smaller dict and find matches in larger dict.
        Write three tsv outputs.
            - differences between matches
            - unique to smaller dict
            - unique to larger dict
        Return three dictionaries for plotting
            - differences between matches
            - unique to smaller dict
            - unique to larger dict
    '''
    # initialize dictionaries
    diffs = {'x': [], 'y': []}
    uFast = {'x': [], 'y': []}
    uBlast = {'x': [], 'y': []}

    # initialize output files
    diffout = f'{outpre}_diffs.tsv'
    ublastout = f'{outpre}_uBlasts.tsv'
    ufastout = f'{outpre}_uFasts.tsv'
    dheader = 'Read pair\tx diff\ty diff\tANI\tpID\tsfrac\tafrac\n'
    bheader = 'Read pair\tpID\tafrac\n'
    fheader = 'Read pair\tANI\tsfrac\n'

    # read through BlastN and find matching pairs in fastANI
    with open(diffout, 'w') as D, open(ublastout, 'w') as B:
        # write headers
        D.write(dheader)
        B.write(bheader)
        for pair, bvalues in BlastN.items():
            if pair in fastANI:
                fvalues = fastANI.pop(pair) # find and remove
                x_diff = fvalues[0] - bvalues[0]
                y_diff = fvalues[1] - bvalues[1]
                lineout = (
                        f'{pair}\t{x_diff}\t{y_diff}\t{fvalues[0]}\t'
                        f'{bvalues[0]}\t{fvalues[1]}\t{bvalues[1]}\n'
                        )
                D.write(lineout)
                diffs['x'].append(x_diff)
                diffs['y'].append(y_diff)
            # if the pair is not in fastANI write it to blastN unique    
            else:
                lineout = f'{pair}\t{bvalues[0]}\t{bvalues[1]}\n'
                B.write(lineout)
                uBlast['x'].append(bvalues[0])
                uBlast['y'].append(bvalues[1])

    # fastANI unique remain. write to file and add to dict
    with open(ufastout, 'w') as F:
        F.write(fheader)
        for pair, fvalues in fastANI.items():
            lineout = f'{pair}\t{fvalues[0]}\t{fvalues[1]}\n'
            F.write(lineout)
            uFast['x'].append(fvalues[0])
            uFast['y'].append(fvalues[1])

    # package the data
    data = {
            'fastANI-blastN': diffs,
            'fastANI_unique': uFast,
            'blastN_unique': uBlast
            }

    return data


def build_hist_plot(xydict, bins, name, title, outpre):

    # Set Colors and markers
    grid_color = '#d9d9d9'
    main_color = '#933b41'
    second_color = '#737373'
    vline_color = '#000000'
    color = '#252525'
    marker = '.' #'o'

    # setup plot
    gg = sns.JointGrid(x=xydict['x'], y=xydict['y'])

    # x margin hist plot
    sns.histplot(
            x=xydict['x'],
            ax=gg.ax_marg_x,
            legend=False,
            color=color,
            stat='probability'
            )
    # y margin hist plot
    sns.histplot(
            y=xydict['y'],
            ax=gg.ax_marg_y,
            legend=False,
            color=color,
            stat='probability'
            )
    # main plot
    gg.ax_joint.plot(
                    xydict['x'],
                    xydict['y'],
                    marker,
                    ms=2,
                    alpha=0.1,
                    color=color,
                    rasterized=True
                    )

    n = len(xydict['x'])
    gg.ax_joint.text(
        0.85, 0.90, f'n = {n}',
        fontsize=8, transform=gg.ax_joint.transAxes
        )

    # plot title, labels, text
    gg.ax_marg_x.set_title(f'{title}: {name}', fontsize=18, y=0.97)
    if name == 'fastANI-blastN':
        gg.ax_joint.set_ylabel('Shared_frac - Aligned_frac', fontsize=12)
        gg.ax_joint.set_xlabel('ANI - pID', fontsize=12)
    elif name == 'fastANI_unique':
        gg.ax_joint.set_ylabel('Shared fraction', fontsize=12)
        gg.ax_joint.set_xlabel('ANI (%)', fontsize=12)
    else:
        gg.ax_joint.set_ylabel('Aligned fraction', fontsize=12)
        gg.ax_joint.set_xlabel('pID (%)', fontsize=12)

    # set the axis parameters / style
    
    xmin = math.floor(min(xydict['x']))
    xmax = math.ceil(max(xydict['x']))
    xstep = (xmax - xmin) / 10
    hstep = xstep/10
    gg.ax_joint.set_xticks(np.arange(xmin, xmax+hstep, xstep))
    gg.ax_joint.set_xlim(left=xmin-hstep, right=xmax+hstep)

    ymin = math.floor(min(xydict['y'])*100)/100
    ymax = math.ceil(max(xydict['y'])*100)/100
    ystep = (ymax - ymin) / 10
    vstep = ystep/10
    gg.ax_joint.set_yticks(np.arange(ymin, ymax+vstep, ystep))
    gg.ax_joint.set_ylim(bottom=ymin-vstep, top=ymax+vstep)

    gg.ax_joint.tick_params(axis='both', labelsize=12)
    gg.ax_joint.tick_params(
        axis='both', which='major', direction='inout', color='k',
        width=2, length=6, bottom=True, left=True, zorder=3
        )
    # set grid style
    gg.ax_joint.yaxis.grid(
        which="major", color='#d9d9d9', linestyle='--', linewidth=1
        )
    gg.ax_joint.xaxis.grid(
        which="major", color='#d9d9d9', linestyle='--', linewidth=1
        )
    gg.ax_joint.set_axisbelow(True)
    gg.ax_joint.add_artist(TickRedrawer())

    # adjust layout, save, and close
    gg.fig.set_figwidth(8)
    gg.fig.set_figheight(6)
    plt.subplots_adjust(left=0.15, right=0.98, top=0.90, bottom=0.10)
    #plt.tight_layout()
    plt.savefig(f'{outpre}_{name}_plot.pdf')
    plt.close()

    return True


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-fa', '--input_fastani_file',
        help='Please specify the input fastANI file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-k', '--fastani_file_key',
        help='Please specify file to read name key for fastANI!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-bn', '--input_blastn_file',
        help='Please specify the input tabular BlastN file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--output_file_prefix',
        help='Please specify the output file prefix!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-t', '--plot_title',
        help='(OPTIONAL) Please specify the plot title!',
        metavar='',
        type=str,
        required=False
        )
    parser.add_argument(
        '-b', '--number_of_bins',
        help='(OPTIONAL) Specify the number of bins (Default = 30).',
        metavar='',
        type=int,
        required=False,
        default=30
        )
    parser.add_argument(
        '-min_ani', '--minimum_ani_value',
        help='(OPTIONAL) Specify the minimum ani value (Default = 70).',
        metavar='',
        type=int,
        required=False,
        default=70
        )
    parser.add_argument(
        '-min_pid', '--minimum_pid_value',
        help='(OPTIONAL) Specify the minimum pID value (Default = 70).',
        metavar='',
        type=int,
        required=False,
        default=70
        )
    parser.add_argument(
        '-min_frac', '--minimum_aligned_frac',
        help='(OPTIONAL) Specify the minimum aligned_frac (Default = 0.5).',
        metavar='',
        type=float,
        required=False,
        default=0.5
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\n\nRunning Script...\n')

    # define parameters
    fastin = args['input_fastani_file']
    key_file = args['fastani_file_key']
    blastin = args['input_blastn_file']
    outpre = args['output_file_prefix']
    title = args['plot_title']
    bins = args['number_of_bins']
    min_ani = args['minimum_ani_value']
    min_pid = args['minimum_pid_value']
    min_frac = args['minimum_aligned_frac']

    if not title:
        fa = fastin.split('/')[-1].split('.')[0]
        bn = blastin.split('/')[-1].split('.')[0]
        title = f'{fa} vs. {bn}'

    # get the read name to file key for fastANI
    print('\nParsing file name to read name key ...\n')
    fast_key = file_to_read_key(key_file)

    # parse the fastANI file - returns dict {read_pair: [ANI, shared_fraction]}
    print('\nParsing fastANI file ...\n')
    fastANI = parse_fastani(fastin, min_ani, min_frac, fast_key)

    # parse the BlastN file - returns dict {read_pair: [pID, aligned_fraction]}
    print('\nParsing the tabular BlastN file ...\n')
    BlastN = parse_blastn(blastin, min_pid, min_frac)

    # compare the two
    data = compute_difference(fastANI, BlastN, outpre)

    # build some plots
    print('\nBuilding plots ...')
    for name, xydict in data.items():
        if len(xydict['x']) > 1:
            _ = build_hist_plot(xydict, bins, name, title, outpre)

    print(f'\n\nComplete success space cadet!! Finished without errors.\n\n')
    
if __name__ == "__main__":
    main()

