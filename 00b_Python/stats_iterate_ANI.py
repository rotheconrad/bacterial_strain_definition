#!/usr/bin/env python3

'''
-------------------------------------------
Author :: Dorian J Feistel
Email :: dfeistel3@gatech.edu
GitHub :: https://github.com/dfesitel3
Date Created :: May 2020
License :: GNU GPLv3
Copyright 2022 Dorian J Feistel
All rights reserved
-------------------------------------------

This script will calculating total data points, precision, recall, F1-score and accuracy and print to stdout
usage:
python3 stats_iterate_ANI.py <Full_ANI_MLST.ani> <ANI_min> > <output>
'''

import time
import sys

def getANI(in_ani):
    '''read the ANI fiel and store are a list (probably not the most memory efficient!)'''
    ST_dict = dict()
    fh = open(in_ani, 'r')
    next(fh) #skip header
    for line in fh:
        line = line.strip().split()
        if line[5] not in ST_dict:
            ST_dict[line[5]] = [line]
        else:
            ST_dict[line[5]].append(line)
    return ST_dict


def get_table(ST, ANI_min):
    '''Determine TN, FN, TP, & TN'''
    table_dict = dict()
    while ANI_min <= 99.99:
        TP = 0
        TN = 0
        FN = 0
        FP = 0
        for ele in ST_dict[ST]:
            if ele[6] == ST and float(ele[2]) >= ANI_min:
                TP += 1
            elif ele[6] != ST and float(ele[2]) < ANI_min:
                TN += 1
            elif ele[6] != ST and float(ele[2]) >= ANI_min:
                FP += 1
            elif ele[6] == ST and float(ele[2]) < ANI_min:
                FN += 1
        table_dict[ST + '_' + str(round(ANI_min, 2))] = [TP, TN, FP, FN]
        ANI_min += 0.1
    return table_dict

def calulate_stats(table_dict, ANI_min):
    '''calulate precision and recall for a given ST'''
    ST_keys = table_dict.keys()
    for ST in ST_keys:
        '''[TP, TN, FP, FN]'''
        confusion_matrix = table_dict[ST]
        TP = confusion_matrix[0]
        TN = confusion_matrix[1]
        FP = confusion_matrix[2]
        FN = confusion_matrix[3]
        data_points = TP + TN + FP + FN
        Precision = round(TP / (TP + FP) * 100, 2)
        Recall = round(TP / (TP + FN) * 100, 2)
        F1_Score = round((2 * Precision * Recall) / (Precision + Recall), 2)
        Accuracy = round((TP + TN) / (TP + TN + FN + FP) * 100, 2)
        print(f'{ST.split("_")[0]}\t{ST.split("_")[1]}\t{data_points}\t{Precision}\t{Recall}\t{F1_Score}\t{Accuracy}')
    return None

if __name__ == "__main__":

    start_time = time.time()
    ANI_infile = sys.argv[1]
    ANI_min = float(sys.argv[2])
    ST_dict = getANI(ANI_infile)
    ST_keys = ST_dict.keys()
    print(f'ST\tANI\tData-Points\tPrecision\tRecall\tF1-Score\tAccuracy')
    for keys in ST_keys:
        internal_time = time.time()
        print(f'Getting statistics for {keys}', file = sys.stderr)
        current_time = round(time.time() - start_time, 2)
        table_dict = get_table(keys, ANI_min)
        calulate_stats(table_dict, ANI_min)
        end_time = round(time.time() - internal_time, 2)
        print(f'time spent: {end_time} sec', file = sys.stderr)
    final_time = round(time.time() - start_time, 2)
    print(f'Fin...\n total time: {final_time} sec', file = sys.stderr)
