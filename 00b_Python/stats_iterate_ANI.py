#!/usr/bin/env python3

'''

this script if for calulating the precision and recall of sequence types (STs) using ANI values
but it will iterate over some staring minimum ANI and up to 100 and give the results

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


def get_table(ST, ANI_min):
    '''Determine TN, FN, TP, & TN'''
    table_dict = dict()
    while ANI_min <= 99.9:
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
        table_dict[ST + '_' + str(round(ANI_min, 1))] = [TP, TN, FP, FN]
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
        try:
            Precision = round(TP / (TP + FP) * 100, 2)
            Recall = round(TP / (TP + FN) * 100, 2)
            F1_Score = round((2 * Precision * Recall) / (Precision + Recall), 2)
            Accuracy = round((TP + TN) / (TP + TN + FN + FP) * 100, 2)
            print(f'{ST.split("_")[0]}\t{ST.split("_")[1]}\t{data_points}\t{Precision}\t{Recall}\t{F1_Score}\t{Accuracy}')
        except ZeroDivisionError:
            continue
    return None

if __name__ == "__main__":

    ANI_infile = sys.argv[1]
    ANI_min = float(sys.argv[2])
    ST_dict = getANI(ANI_infile)
    ST_keys = ST_dict.keys()
    print(f'ST\tANI\tData-Points\tPrecision\tRecall\tF1-Score\tAccuracy')
    for keys in ST_keys:
        table_dict = get_table(keys, ANI_min)
        calulate_stats(table_dict, ANI_min)
