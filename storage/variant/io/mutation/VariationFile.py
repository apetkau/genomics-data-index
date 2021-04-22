from __future__ import annotations

import uuid
import logging
import tempfile
from pathlib import Path
from typing import List, Tuple
import multiprocessing as mp

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from storage.variant.util import execute_commands

logger = logging.getLogger(__name__)


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
    def union_all_files(cls, variant_files: List[Path], include_expression: str = 'TYPE="SNP"',
                        ncores=2) -> pd.DataFrame:

        if len(variant_files) == 0:
            raise Exception('Cannot take union of 0 files')
        elif len(variant_files) == 1:
            return cls.read_single_file(variant_files[0], include_expression=include_expression)
        else:
            return cls._union_multiple_files(variant_files, include_expression=include_expression,
                                             ncores=ncores)

    @classmethod
    def read_single_file(cls, variant_file: Path, include_expression: str = 'TYPE="SNP"') -> pd.DataFrame:
        with tempfile.TemporaryDirectory() as tmp_dir:
            output_file = Path(tmp_dir) / 'output.tsv'

            command = ['bcftools', 'query', '-f', '%CHROM\t%POS\t%REF\t%ALT\t1\n']
            if include_expression is not None:
                command.extend(['-i', include_expression])
            command.extend(['-o', str(output_file), str(variant_file)])

            execute_commands([
                command
            ])
            var_df = pd.read_csv(output_file, sep='\t', dtype=str,
                                 names=['CHROM', 'POS', 'REF', 'ALT', 'INDEXES'])
            var_df['POS'] = var_df['POS'].astype(int)

            # Should only be one entry in a single file, so count is always 1
            var_df['COUNT'] = 1
            var_df = var_df.drop(columns='INDEXES')

            return var_df.sort_values(['CHROM', 'POS'])

    @classmethod
    def _union_multiple_files(cls, variant_files: List[Path], include_expression: str = 'TYPE="SNP"',
                              ncores: int = 1, batch_size: int = 50) -> pd.DataFrame:
        with tempfile.TemporaryDirectory() as tmp_dir_str:
            # First batch up files to be processed by bcftools isec (cannot handle very long lists of files)
            tmp_dir = Path(tmp_dir_str)
            batches = []
            tmp_file = Path(tmp_dir, str(uuid.uuid1()))
            batch = BcfToolsUnionExecutor(nthreads=1, include_expression=include_expression,
                                          out_file=tmp_file)
            batches.append(batch)
            num_in_batch = 0
            for file in variant_files:
                batch.add(file)
                num_in_batch += 1
                # start new batch
                if num_in_batch == batch_size:
                    tmp_file = Path(tmp_dir, str(uuid.uuid1()))
                    batch = BcfToolsUnionExecutor(nthreads=1, include_expression=include_expression,
                                                  out_file=tmp_file)
                    batches.append(batch)
                    num_in_batch = 0

            # Execute batches in a process pool
            logger.debug(f'Started processing {len(variant_files)} VCF/BCF files in batches of {batch_size}')
            with mp.Pool(ncores) as pool:
                union_df_results = []

                union_dfs = pool.map(cls._single_batch_union, batches)

                for union_df in union_dfs:
                    union_df_results.append(union_df)

                var_df = pd.concat(union_df_results).groupby('ID').agg({
                    'ID': 'first',
                    'CHROM': 'first',
                    'POS': 'first',
                    'REF': 'first',
                    'ALT': 'first',
                    'COUNT': 'sum',
                })
                var_df['POS'] = var_df['POS'].astype(int)

                logger.debug(f'Finished processing {len(variant_files)} VCF/BCF files')
                return var_df.sort_values(['CHROM', 'POS'])

    @classmethod
    def _single_batch_union(cls, union_executor: BcfToolsUnionExecutor) -> pd.DataFrame:
        return union_executor.union()


class BcfToolsUnionExecutor:

    def __init__(self, nthreads: int, include_expression: str, out_file: Path):
        self._files = []
        self._nthreads = nthreads
        self._include_expression = include_expression
        self._out_file = out_file

    def add(self, file: Path):
        self._files.append(file)

    def union(self) -> pd.DataFrame:
        command = ['bcftools', 'isec']
        if self._include_expression is not None:
            command.extend(['-i', self._include_expression])
        command.extend(['-c', 'none', '-n', '+1', '--threads', f'{self._nthreads}', '-o', str(self._out_file)])
        for file in self._files:
            command.append(str(file))

        execute_commands([
            command
        ])

        var_df = pd.read_csv(self._out_file, sep='\t', dtype=str,
                             names=['CHROM', 'POS', 'REF', 'ALT', 'INDEXES'])
        # Count occurences of '1' character in INDEXES which represents number of samples with this mutation
        var_df['COUNT'] = var_df['INDEXES'].apply(lambda x: x.count('1'))
        var_df = var_df.drop(columns=['INDEXES'])

        # Added this ID so I can group by later
        var_df['ID'] = var_df.apply(lambda x: f'{x["CHROM"]}:{x["POS"]}:{x["REF"]}:{x["ALT"]}', axis='columns')

        return var_df
