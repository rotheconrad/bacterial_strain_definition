#!/usr/bin/env python3

import sys

def open_MLST(inMLST):
    mlst_dict = {}
    with open(inMLST, 'r') as fh:
        for line in fh:
            line = line.strip().split()
            mlst_dict[line[0]] = line[2]
    return mlst_dict

def open_ANI_file(inANI, mlst_dict):
    with open(inANI, 'r') as fh, open(inANI + "_stats_file.ani", "w") as wf:
        wf.write(f"query\treference\tANI\tfragments\ttotal\tquery_ST\tST\tproportion\n")
        for line in fh:
            line = line.strip().split()
            if line[0] == line[1]: continue
            prop = round(int(line[3]) / int(line[4]), 4)
            output = '\t'.join(line)
            wf.write(f'{output}\tST-{mlst_dict[line[0]]}\tST-{mlst_dict[line[1]]}\t{prop}\n')
    return None

if __name__ == "__main__":

    MLST_dict = open_MLST(sys.argv[1])
    open_ANI_file(sys.argv[2], MLST_dict)



