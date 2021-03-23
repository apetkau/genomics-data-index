from typing import List

import tempfile
from pathlib import Path

from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from storage.variant.io import execute_commands


class VariationFile:

    def __init__(self, file: Path):
        self._file = file

    def write(self, output: Path, file_type: str = 'bcf') -> Path:
        if file_type == 'bcf':
            execute_commands([
                ['bcftools', 'view', str(self._file), '-o', str(output), '-O', 'b'],
                ['bcftools', 'index', str(output)]
            ])
            return output
        else:
            raise Exception(f'Invalid file_type=[{file_type}]')

    def consensus(self, reference_file: Path, mask_file: Path, include_expression='TYPE="snp"',
                  mask_with: str = 'N') -> List[SeqRecord]:
        with tempfile.NamedTemporaryFile() as out_f:
            execute_commands([
                ['bcftools', 'consensus', '--fasta-ref', str(reference_file),
                 '--mask', str(mask_file), '--mask-with', mask_with,
                 '--include', include_expression,
                 '--output', str(out_f.name),
                 str(self._file)]
            ])
            sequences = list(SeqIO.parse(out_f.name, 'fasta'))
            return sequences
