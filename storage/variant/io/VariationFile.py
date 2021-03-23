from typing import List

import tempfile
from pathlib import Path
import subprocess

from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


class VariationFile:

    def __init__(self, file: Path):
        self._file = file

    def execute_commands(self, commands: List[List[str]]):
        try:
            for command in commands:
                subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, text=True)
        except subprocess.CalledProcessError as e:
            err_msg = str(e.stderr.strip())
            raise Exception(f'Could not run [{" ".join(e.cmd)}] on original file [{self._file}]: error {err_msg}')

    def write(self, output: Path, file_type: str = 'bcf') -> Path:
        if file_type == 'bcf':
            self.execute_commands([
                ['bcftools', 'view', str(self._file), '-o', str(output), '-O', 'b'],
                ['bcftools', 'index', str(output)]
            ])
            return output
        else:
            raise Exception(f'Invalid file_type=[{file_type}]')

    def consensus(self, reference_file: Path, mask_file: Path, include_expression='TYPE="snp"') -> List[SeqRecord]:
        with tempfile.NamedTemporaryFile() as out_f:
            self.execute_commands([
                ['bcftools', 'consensus', '--fasta-ref', str(reference_file),
                 '--mask', str(mask_file), '--mask-with', 'N',
                 '--include', include_expression,
                 '--output', str(out_f.name),
                 str(self._file)]
            ])
            sequences = list(SeqIO.parse(out_f.name, 'fasta'))
            return sequences
