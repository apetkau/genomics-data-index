from pathlib import Path

from genomics_data_index.storage.io.mutation.SequenceFile import SequenceFile
from genomics_data_index.storage.MaskedGenomicRegions import MaskedGenomicRegions


# Generates a mask bed file from a snippy output directory


snippy_dir = Path(snakemake.input.snippy_dir)
output_mask = Path(snakemake.output.mask)

consensus_file = snippy_dir / 'snps.aligned.fa'

name, consensus_records = SequenceFile(consensus_file).parse_sequence_file()
mask = MaskedGenomicRegions.from_sequences(consensus_records)

mask.write(output_mask)
