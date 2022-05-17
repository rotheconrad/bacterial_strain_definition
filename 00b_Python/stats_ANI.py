#!/usr/bin/env python3

'''

this script if for calculating the precision and recall of sequence types (STs) using ANI values
it will total data points, precision, recall, F1-score and accuracy

stats_ANI.py <parsed_ANI_file> <ANI_value> > <output_file>

Note: output is printed to stdout!

'''

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
            if ele[6] == ST and float(ele[2]) > ANI_min:
                TP += 1
            elif ele[6] != ST and float(ele[2]) < ANI_min:
                TN += 1
            elif ele[6] != ST and float(ele[2]) > ANI_min:
                FP += 1
            elif ele[6] == ST and float(ele[2]) < ANI_min:
                FN += 1
        table_dict[ST] = [TP, TN, FP, FN]
    return table_dict

def calulate_stats(table_dict, ANI_min):
    '''calulate precision and recall for a given ST'''
    print(f'Precision and recall for a given ST at a minimum ANI of {ANI_min}')
    print(f'ST\tData-points\tPrecision\tRecall\tF1-Score\tAccuracy')
    ST_keys = table_dict.keys()
    for ST in ST_keys:
        '''[TP, TN, FP, FN]'''
        confusion_matrix = table_dict[ST]
        TP = confusion_matrix[0]
        TN = confusion_matrix[1]
        FP = confusion_matrix[2]
        FN = confusion_matrix[3]
        total_data = TP + TN + FP + FN
        Precision = round(TP / (TP + FP) * 100, 2)
        Recall = round(TP / (TP + FN) * 100, 2)
        F1_Score = round((2 * Precision * Recall) / (Precision + Recall), 2)
        Accuracy = round((TP + TN) / (TP + TN + FN + FP) * 100, 2)
        print(f'{ST}\t{total_data}\t{Precision}\t{Recall}\t{F1_Score}\t{Accuracy}')
    return None

if __name__ == "__main__":

    ANI_infile = sys.argv[1]
    ANI_min = float(sys.argv[2])
    ST_dict = getANI(ANI_infile)
    table_dict = get_table(ST_dict, ANI_min)
    calulate_stats(table_dict, ANI_min)
