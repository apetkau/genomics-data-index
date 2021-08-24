#!/usr/bin/env python
# Prints the length of all intervals in a bed file
# Usage: cat file.bed | print-length-bed.py 

import sys
from pathlib import Path

from pybedtools import BedTool

lines = ''.join(sys.stdin.readlines())
bed = BedTool(lines, from_string=True)

total = 0
for i in bed:
    total += len(i)

print(total)

