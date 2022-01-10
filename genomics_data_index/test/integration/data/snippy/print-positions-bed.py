#!/usr/bin/env python
# Prints the sequential positions of all intervals in a bed file (1-based coordinates)
# Usage: cat file.bed | print-positions-bed.py 

import sys

from pybedtools import BedTool

lines = ''.join(sys.stdin.readlines())
bed = BedTool(lines, from_string=True)

for i in bed:
    start = i.start + 1
    stop = i.stop + 1
    for pos in range(start, stop):
        print(f'{i.chrom}\t{pos}')
