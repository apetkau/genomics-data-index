#!/usr/bin/env python
# Converts snippy aligned.fa and VCF to a BED mask file

import sys
from pathlib import Path

from pybedtools import BedTool

from genomics_data_index.storage.MaskedGenomicRegions import MaskedGenomicRegions
from genomics_data_index.storage.io.mutation.SequenceFile import SequenceFile

if len(sys.argv) != 3:
    print(f"Usage: {sys.argv[0]} [snps.aligned.fa] [snps.vcf.gz]\n")
    sys.exit(1)

snps_aligned = Path(sys.argv[1])
snps_vcf = Path(sys.argv[2])

sf = SequenceFile(snps_aligned)
vcf_bed = BedTool(snps_vcf).merge()
name, records = list(sf.parse_sequence_file())
regions = MaskedGenomicRegions.from_sequences(records)
regions_minus_vcf = MaskedGenomicRegions(regions._mask.subtract(vcf_bed))

output_file = Path(name + '.bed.gz')
regions.write(output_file)
