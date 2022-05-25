#!/usr/bin/env python3

'''
Use this script for parsing and preparing FastANI data for plotting Proportion of Fragments vs ANI

usage:
python3 parse_fastANI_data.py <ANI_file> <MLST_file>

-------------------------------------------
Author :: Dorian J Feistel
Email :: dfeistel3@gatech.edu
GitHub :: https://github.com/dfesitel3
Date Created :: May 2020
License :: GNU GPLv3
Copyright 2022 Dorian J Feistel
All rights reserved
-------------------------------------------
'''

import sys
import time

def open_MLST_file(in_mlst):
    mlst_dict = dict()
    fh = open(in_mlst, 'r')
    for line in fh:
        line = line.strip().split()
        genome = line[0]
        ST_num = line[2]
        ST_str = 'ST-' + ST_num
        mlst_dict[genome] = ST_str
    fh.close()
    return mlst_dict

def open_ANI_file(infile, mlst_dict):
    '''open ANI file, removed self-matches and smaller of two genomes,
    calculate proportion of fragments, create Dictionary for writing'''
    print('Removing self-matches, keeping larger of two genomes, and calculating proportion of fragments')
    ANI_dict = dict()
    fh2 = open(infile, 'r')
    count = 0
    for line in fh2:
        line = line.strip().split()
        query = line[0]
        ref = line[1]
        frags = int(line[3])
        total_frags = int(line[4])

        if query != ref:
            '''add mlst for query and reference'''
            line.append(mlst_dict[query])
            line.append(mlst_dict[ref])
            '''calculate fragment proportion'''
            prop = round(frags/total_frags, 4)
            line.append(str(prop))
            '''removed smaller genome'''
            names = [query, ref]
            names.sort()
            gname = '-'.join(names)
            if gname not in ANI_dict:
                ANI_dict[gname] = line
            else:
                if ANI_dict[gname][4] < line[4]:
                    ANI_dict[gname] = line
            count += 1
            if count % 1000000 == 0:
                print(f'Processed {count} lines from {infile}')
    fh2.close()
    return ANI_dict

def write_output(ANI_dict, ANI_raw):
    outname = ANI_raw + '_prepaired.ani'
    print(f'Writing parsed data to {outname}')
    wf = open(outname, 'w')
    wf.write(f'query\treference\tANI\tfragments\ttotal\tquery_ST\tST\tproportion\n')
    for k, v in ANI_dict.items():
        wf.write('%s\n' % ('\t'.join(v)))
    wf.close()
    return print(f'Finished parsing!')

def main():
    t0 = time.time()
    ANI_raw = sys.argv[1]
    mlst_raw = sys.argv[2]
    mlst_dict = open_MLST_file(mlst_raw)
    ANI_dict = open_ANI_file(ANI_raw, mlst_dict)
    write_output(ANI_dict, ANI_raw)
    tf = time.time() - t0
    print(f'Processed in {tf}')
    return None

if __name__ == "__main__":
    main()
