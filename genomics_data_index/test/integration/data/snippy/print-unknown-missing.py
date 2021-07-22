#!/usr/bin/env python
# Used to print out unknown/missing ('N' or '-') in snippy alignment files for debugging

import sys

from Bio import SeqIO

if len(sys.argv) != 2:
    print(f"Usage: {sys.argv[0]} [snps.aligned.fa]\n")
    sys.exit(1)

with open(sys.argv[1]) as handle:
    for record in SeqIO.parse(handle, "fasta"):
        name = record.id
        for i, c in enumerate(record.seq):
            c = c.upper()
            if c == "N" or c == "-":
                pos = i + 1
                print(pos)
