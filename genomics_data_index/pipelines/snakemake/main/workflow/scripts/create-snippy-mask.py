from pathlib import Path

from genomics_data_index.storage.MaskedGenomicRegions import MaskedGenomicRegions
from genomics_data_index.storage.io.mutation.SequenceFile import SequenceFile

# Generates a mask bed file from a snippy output directory


snippy_dir = Path(snakemake.input.snippy_dir)
output_mask = Path(snakemake.output.mask)

consensus_file = snippy_dir / 'snps.aligned.fa'
vcf_file = snippy_dir / 'snps.vcf.gz'

name, consensus_records = SequenceFile(consensus_file).parse_sequence_file()
mask = MaskedGenomicRegions.from_sequences(consensus_records)
vcf_regions = MaskedGenomicRegions.from_vcf_file(vcf_file)
mask = mask.subtract(vcf_regions)

mask.write(output_mask)
