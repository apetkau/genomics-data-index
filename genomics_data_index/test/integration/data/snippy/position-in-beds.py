#!/usr/bin/env python
# Prints if the given position (1-based, like in VCF files) is in the given bed files

import sys
from pathlib import Path

from genomics_data_index.storage.MaskedGenomicRegions import MaskedGenomicRegions

if len(sys.argv) < 2:
    print(f"Usage: {sys.argv[0]} [position] [file.bed] ...\n")
    sys.exit(1)

seq_position = sys.argv[1]
if ':' in seq_position:
    sequence, position = seq_position.split(':')
    position = int(position)
else:
    sequence = 'reference'
    position = int(seq_position)

for file in sys.argv[2:]:
    mask = MaskedGenomicRegions.from_file(Path(file))
    if mask.contains(sequence=sequence, position=position, start_position_index='1'):
        print(f'{file} contains {sequence}:{position}')
#    else:
#        print(f'{file} does not contain {sequence}:{position}')
