#!/usr/bin/env python3

import sys

data_dict = {}
with open(sys.argv[1], 'r') as fh:
	for line in fh:
		line = line.strip().split()
		data_dict[line[2]] = data_dict.get(line[2], 0) + 1


for i in sorted(data_dict.keys()):
	if data_dict[i] >= 10:
		print(f'{i}\t{data_dict[i]}')
