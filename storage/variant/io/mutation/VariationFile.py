import tempfile
from pathlib import Path
from typing import List, Tuple

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from storage.variant.util import execute_commands


class VariationFile:

    def __init__(self, file: Path):
        self._file = file

    def write(self, output: Path) -> Tuple[Path, Path]:
        if output.suffix == '.bcf':
            output_type = 'b'
        elif str(output).endswith('.vcf.gz'):
            output_type = 'z'
        else:
            raise Exception((f'Invalid file type [{output.suffix}] for output=[{output}]. '
                             'Must be one of [".bcf", ".vcf.bz"]'))

        output_index = Path(str(output) + '.csi')

        execute_commands([
            ['bcftools', 'plugin', 'fill-tags', str(self._file), '-O', output_type, '-o', str(output),
             '--', '-t', 'TYPE'],
            ['bcftools', 'index', str(output)]
        ])
        return output, output_index

    def consensus(self, reference_file: Path, mask_file: Path = None, include_expression='TYPE="SNP"',
                  mask_with: str = 'N') -> List[SeqRecord]:
        with tempfile.NamedTemporaryFile() as out_f:
            command = ['bcftools', 'consensus', '--fasta-ref', str(reference_file)]
            if mask_file is not None:
                command.extend(['--mask', str(mask_file), '--mask-with', mask_with])
            command.extend([
                '--include', include_expression,
                '--output', str(out_f.name),
                str(self._file)
            ])
            execute_commands([command])
            sequences = list(SeqIO.parse(out_f.name, 'fasta'))
            return sequences

    @classmethod
    def union_all_files(cls, variant_files: List[Path], include_expression: str = 'TYPE="SNP"') -> pd.DataFrame:
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
                                 names=['CHROM', 'POS', 'REF', 'ALT', 'INDEXES'])
            var_df['POS'] = var_df['POS'].astype(int)
            return var_df.sort_values(['CHROM', 'POS'])
