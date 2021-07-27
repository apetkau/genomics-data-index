from __future__ import annotations

import logging
import multiprocessing as mp
import os
import tempfile
import uuid
from pathlib import Path
from typing import List, Tuple, Optional

import pandas as pd
import vcf
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from genomics_data_index.storage.io.mutation.SnpEffDatabase import SnpEffDatabase
from genomics_data_index.storage.io.mutation.VcfSnpEffAnnotationParser import VcfSnpEffAnnotationParser
from genomics_data_index.storage.model.NucleotideMutationTranslater import NucleotideMutationTranslater
from genomics_data_index.storage.util import execute_commands

logger = logging.getLogger(__name__)


def translate_to_mutation_id(x: pd.Series) -> str:
    return NucleotideMutationTranslater.to_spdi(sequence_name=x['CHROM'],
                                                position=x['POS'],
                                                ref=x['REF'],
                                                alt=x['ALT'],
                                                convert_deletion=False)


class VariationFile:
    VCF_FRAME_COLUMNS = ['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'FILE', 'VARIANT_ID']

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

    def _fix_df_columns(self, vcf_df: pd.DataFrame) -> pd.DataFrame:
        # If no data, I still want certain column names so that rest of code still works
        if len(vcf_df) == 0:
            vcf_df = pd.DataFrame(columns=['CHROM', 'POS', 'REF', 'ALT', 'INFO'])

        return vcf_df.loc[:, ['CHROM', 'POS', 'REF', 'ALT', 'INFO']]

    def _fix_alt(self, element: List[str]) -> str:
        """
        Fix up the alternative string as the pyVCF package does not return them as a string.
        :param element: The element to fix.
        :return: The fixed element.
        """
        return str(element[0])

    def _fix_ref(self, element: List[str]) -> str:
        """
        Fix up the reference string as the pyVCF package does not return them as a string.
        :param element: The element to fix.
        :return: The fixed element.
        """
        return str(element)

    def _get_type(self, vcf_df: pd.DataFrame) -> pd.Series:
        return vcf_df['INFO'].map(lambda x: x['TYPE'][0])

    def _drop_extra_columns(self, vcf_df: pd.DataFrame) -> pd.DataFrame:
        return vcf_df.drop('INFO', axis='columns')

    def read_features(self, sample_name: str, snpeff_parser: VcfSnpEffAnnotationParser) -> pd.DataFrame:
        reader = vcf.Reader(filename=str(self.file))
        df = pd.DataFrame([vars(r) for r in reader])
        out = self._fix_df_columns(df)

        out['ALT'] = out['ALT'].map(self._fix_alt)
        out['REF'] = out['REF'].map(self._fix_ref)
        out['TYPE'] = self._get_type(out)

        snpeff_headers = snpeff_parser.parse_annotation_headers(vcf_info=reader.infos)
        ann_df = snpeff_parser.parse_annotation_entries(vcf_ann_headers=snpeff_headers, vcf_df=out)
        out = out.merge(ann_df, how='left', left_index=True, right_index=True)

        out = self._drop_extra_columns(out)

        out['FILE'] = os.path.basename(self.file)
        cols = out.columns.tolist()
        out['SAMPLE'] = sample_name
        out = out.reindex(columns=['SAMPLE'] + cols)
        out = snpeff_parser.select_variant_annotations(out)
        return out.loc[:, self.VCF_FRAME_COLUMNS + snpeff_parser.ANNOTATION_COLUMNS]

    @property
    def file(self) -> Path:
        return self._file

    def annotate(self, snpeff_database: SnpEffDatabase, annotated_vcf: Path) -> VariationFile:
        return VariationFile(snpeff_database.annotate(self._file, output_vcf_file=annotated_vcf))

    def consensus(self, reference_file: Path, mask_file: Path = None, include_expression: str = 'TYPE="SNP"',
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
    def union_all_files(cls, variant_files: List[Path], include_expression: Optional[str] = None,
                        ncores=1, batch_size=50) -> pd.DataFrame:
        if len(variant_files) == 0:
            empty_df = pd.DataFrame([], columns=['ID', 'CHROM', 'POS', 'REF', 'ALT', 'COUNT'])
            empty_df['POS'] = empty_df['POS'].astype(int)
            return empty_df
        else:
            with tempfile.TemporaryDirectory() as tmp_dir_str:
                # First batch up files to be processed (tools cannot handle too many command-line arguments and
                # bcftools isec does not support reading list of file names from a file)
                # Plus, batch processing means I can distribute over multiple processes
                tmp_dir = Path(tmp_dir_str)
                batches = []
                tmp_file = Path(tmp_dir, str(uuid.uuid1()))
                batch = BcfToolsUnionExecutor(nthreads=1, include_expression=include_expression,
                                              out_file=tmp_file)
                batches.append(batch)
                num_in_batch = 0
                for file in variant_files:
                    batch.add_file(file)
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

    def add_file(self, file: Path):
        self._files.append(file)

    def union(self) -> pd.DataFrame:
        if len(self._files) == 0:
            empty_df = pd.DataFrame([], columns=['ID', 'CHROM', 'POS', 'REF', 'ALT', 'COUNT'])
            empty_df['POS'] = empty_df['POS'].astype(int)
            return empty_df
        elif len(self._files) == 1:
            return self._union_single(self._files[0])
        else:
            return self._union_multiple(self._files)

    def _union_single(self, file: Path) -> pd.Dataframe:
        command = ['bcftools', 'query', '-f', '%CHROM\t%POS\t%REF\t%ALT\n']
        if self._include_expression is not None:
            command.extend(['-i', self._include_expression])
        command.extend(['-o', str(self._out_file), str(file)])

        execute_commands([
            command
        ])
        var_df = pd.read_csv(self._out_file, sep='\t', dtype=str,
                             names=['CHROM', 'POS', 'REF', 'ALT'])
        var_df['POS'] = var_df['POS'].astype(int)

        # Only single file so count is 1
        var_df['COUNT'] = 1

        # Added this ID so I can group by later
        if len(var_df) == 0:
            var_df['ID'] = None
        else:
            var_df['ID'] = var_df.apply(translate_to_mutation_id, axis='columns')

        return var_df.sort_values(['CHROM', 'POS'])

    def _union_multiple(self, files: List[Path]) -> pd.DataFrame:
        command = ['bcftools', 'isec']
        if self._include_expression is not None:
            command.extend(['-i', self._include_expression])
        command.extend(['-c', 'none', '-n', '+1', '--threads', f'{self._nthreads}', '-o', str(self._out_file)])
        for file in files:
            command.append(str(file))

        execute_commands([
            command
        ])

        var_df = pd.read_csv(self._out_file, sep='\t', dtype=str,
                             names=['CHROM', 'POS', 'REF', 'ALT', 'INDEXES'])
        # Count occurences of '1' character in INDEXES which represents number of samples with this mutation
        var_df['COUNT'] = var_df['INDEXES'].apply(lambda x: x.count('1'))
        var_df = var_df.drop(columns=['INDEXES'])

        var_df['POS'] = var_df['POS'].astype(int)

        # Added this ID so I can group by later
        var_df['ID'] = var_df.apply(translate_to_mutation_id, axis='columns')

        return var_df
