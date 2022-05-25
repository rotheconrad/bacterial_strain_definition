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
python3 stats_ANI.py <Full_ANI.ani> <ANI_min_value>
'''
import sys

def getANI(in_ani):
    '''read full ANI file and store as a dict'''
    ST_dict = dict()
    fh = open(in_ani, 'r')
    next(fh) #skip header
    '''create dict with query ST as key and line as values'''
    for line in fh:
        line = line.strip().split()
        ANI = line[2]
        query = line[5]
        reference = line[6]
        if query not in ST_dict:
            ST_dict[query] = [[reference, ANI]]
        else:
            ST_dict[query].append([reference, ANI])
    return ST_dict

def get_table(ST_dict, ANI_min):
    '''Determine TN, FN, TP, & TN'''
    ST_keys = ST_dict.keys()
    table_dict = dict()
    for ST in ST_keys:
        TP = 0
        TN = 0
        FN = 0
        FP = 0
        for ele in ST_dict[ST]:
            ST_compare = ele[0]
            ANI = ele[1]
            if ST_compare == ST and float(ANI) >= ANI_min:
                TP += 1
            elif ST_compare != ST and float(ANI) < ANI_min:
                TN += 1
            elif ST_compare != ST and float(ANI) >= ANI_min:
                FP += 1
            elif ST_compare == ST and float(ANI) < ANI_min:
                FN += 1
        table_dict[ST] = [TP, TN, FP, FN]
    return table_dict

def calulate_stats(table_dict, ANI_min, ST_dict):
    '''calulate precision and recall for a given ST'''
    #print(f'Precision and recall for a given ST at a minimum ANI of {ANI_min}')
    #print(f'ST\tData-points\tPrecision\tRecall\tF1-Score\tAccuracy')
    ST_keys = table_dict.keys()
    for ST in ST_keys:
        '''[TP, TN, FP, FN]'''
        confusion_matrix = table_dict[ST]
        TP = confusion_matrix[0]
        TN = confusion_matrix[1]
        FP = confusion_matrix[2]
        FN = confusion_matrix[3]
        total_data = len(ST_dict[ST])
        Precision = round(TP / (TP + FP) * 100, 2)
        Recall = round(TP / (TP + FN) * 100, 2)
        F1_Score = round((2 * Precision * Recall) / (Precision + Recall), 2)
        Accuracy = round((TP + TN) / (TP + TN + FN + FP) * 100, 2)
        print(f'{ST}\t{total_data}\t{Precision}\t{Recall}\t{F1_Score}\t{Accuracy}')
        #print(f'{TP}\t{TN}\t{FP}\t{FN}')
    return None

if __name__ == "__main__":

    ANI_infile = sys.argv[1]
    ANI_min = float(sys.argv[2])
    ST_dict = getANI(ANI_infile)
    table_dict = get_table(ST_dict, ANI_min)
    calulate_stats(table_dict, ANI_min, ST_dict)
