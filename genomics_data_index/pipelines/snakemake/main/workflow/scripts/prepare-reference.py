from pathlib import Path

from genomics_data_index.storage.io.mutation.SequenceFile import SequenceFile

# Reads and prepares a reference genome by writing it to a FASTA file.

input_reference = Path(snakemake.input.reference)
output_reference = Path(snakemake.output.reference)

sequence_file = SequenceFile(input_reference)
sequence_file.write(output_reference, 'fasta')
