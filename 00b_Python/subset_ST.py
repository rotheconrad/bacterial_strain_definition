#!/usr/bin/env python3

'''
run this script to parse STs after using parse_fastANI_data.py

get_proper_STs.py <parsed_ANI_file> <#,#,#,#>

argument $2 can be one or many digits corresponding to STs in the parsed ANI file, sepereted by a comma e.g. 131,10,11

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

def list_mlst(in_mlst):
	output = ['ST-' + i for i in in_mlst.split(',')]
	return output

def open_paresed_ANI(inANI, mlst_list):
	with open(inANI, 'r') as fh, open(inANI + '_subset_MLSTs.ani', 'w') as wf:
		for line in fh:
			if not line.startswith('query'):
				line = line.strip().split()
				if line[5] in mlst_list and line[6] in mlst_list:
					wf.write('%s\n' % ('\t'.join(line)))
			else:
				line = line.strip().split()
				wf.write('%s\n' % ('\t'.join(line)))
	return None

def main():
	ANI_infile = sys.argv[1]
	MLSTs = sys.argv[2]
	print_mlst = ', '.join(list_mlst(MLSTs))
	print(f'\nParsing {ANI_infile} \n\tLooking for\n\t\t{print_mlst}')
	open_paresed_ANI(ANI_infile, list_mlst(MLSTs))
	print('\nfin...\n')
	return None

if __name__ == "__main__":
	main()


