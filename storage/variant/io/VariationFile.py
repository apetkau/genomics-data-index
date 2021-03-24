from typing import List

import tempfile
from pathlib import Path

import pandas as pd
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

    @classmethod
    def union_all_files(cls, variant_files: List[Path], include_expression: str = 'TYPE="snp"') -> pd.DataFrame:
        with tempfile.TemporaryDirectory() as tmp_dir:
            union_file = Path(tmp_dir) / 'union.tsv'

            if len(variant_files) == 0:
                raise Exception('Cannot take union of 0 files')
            elif len(variant_files) == 1:
                command = ['bcftools', 'query', '-f', '%CHROM\t%POS\t%REF\t%ALT\t1\n', '-i', include_expression,
                           '-o', str(union_file), str(variant_files[0])]
            else:
                command = ['bcftools', 'isec', '-i', include_expression, '-c', 'none',
                           '-n', '+1', '--threads', '2', '-o', str(union_file)]
                for file in variant_files:
                    command.append(str(file))

            execute_commands([
                command
            ])
            var_df = pd.read_csv(union_file, sep='\t', dtype=str,
                               names=['CHROM', 'POS', 'REF', 'ALT', 'INDEXES']).sort_values(['CHROM','POS'])
            var_df['POS'] = var_df['POS'].astype(int)
            return var_df
