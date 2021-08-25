#!/usr/bin/env python
# Converts snippy aligned.fa (and optionally a VCF) to a BED mask file

import sys
from pathlib import Path

from pybedtools import BedTool

from genomics_data_index.storage.MaskedGenomicRegions import MaskedGenomicRegions
from genomics_data_index.storage.io.mutation.SequenceFile import SequenceFile

if len(sys.argv) == 2:
    snps_aligned = Path(sys.argv[1])
    snps_vcf = None
elif len(sys.argv) == 3:
    snps_aligned = Path(sys.argv[1])
    snps_vcf = Path(sys.argv[2])
else:
    print(f"Usage: {sys.argv[0]} snps.aligned.fa [snps.vcf.gz]\n")
    sys.exit(1)

sf = SequenceFile(snps_aligned)
name, records = list(sf.parse_sequence_file())
regions = MaskedGenomicRegions.from_sequences(records)

if snps_vcf is not None:
    name = name + '.minus-vcf'
    vcf_bed = BedTool(snps_vcf).merge()
    regions = MaskedGenomicRegions(regions.mask.subtract(vcf_bed))

output_file = Path(name + '.bed.gz')
regions.write(output_file)
