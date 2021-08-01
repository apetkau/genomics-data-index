#!/usr/bin/env python
# Used to print out the nucleotide and position in the reference genome to fill in information for tests.

import gzip
import sys
from functools import partial

from Bio import SeqIO

if len(sys.argv) != 2:
    print(f"Usage: {sys.argv[0]} [reference.fasta]\n")
    sys.exit(1)
file = sys.argv[1]

_open = partial(gzip.open, mode='rt') if file.endswith('.gz') else open

with _open(file) as handle:
    for record in SeqIO.parse(handle, "fasta"):
        name = record.id
        for i, c in enumerate(record.seq):
            c = c.upper()
            pos = i + 1
            print(f'{pos}\t{c}')
