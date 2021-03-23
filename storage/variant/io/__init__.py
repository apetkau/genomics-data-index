import abc
from typing import Dict, List
from pathlib import Path
import subprocess

import pandas as pd

from storage.variant.MaskedGenomicRegions import MaskedGenomicRegions


def check_variants_table_columns(df: pd.DataFrame) -> None:
    expected_columns = {'SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE'}
    actual_columns = set(df.columns.tolist())
    if not expected_columns.issubset(actual_columns):
        raise Exception('Variants table does not contain expected set of columns. '
                        f'Expected {expected_columns}, actual {actual_columns}')


def execute_commands(commands: List[List[str]]):
    try:
        for command in commands:
            subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, text=True)
    except subprocess.CalledProcessError as e:
        err_msg = str(e.stderr.strip())
        raise Exception(f'Could not run [{" ".join(e.cmd)}]: error {err_msg}')


class VariantsReader(abc.ABC):

    def __init(self):
        pass

    def get_variants_table(self) -> pd.DataFrame:
        variants_df = self._read_variants_table()
        check_variants_table_columns(variants_df)
        return variants_df

    @abc.abstractmethod
    def sample_variant_files(self) -> Dict[str, Path]:
        """
        Gets a dictionary of sample names to variant files to be read by this reader.
        :return: A dictionary of sample names to variant files ('name' => 'file')
        """
        pass

    @abc.abstractmethod
    def samples_list(self) -> List[str]:
        """
        Gets a list of sample names that will be read by this reader.
        :return: A list of sample names that will be read by this reader.
        """
        pass

    @abc.abstractmethod
    def _read_variants_table(self) -> pd.DataFrame:
        pass

    def get_genomic_masked_regions(self) -> Dict[str, MaskedGenomicRegions]:
        return self._read_genomic_masked_regions()

    @abc.abstractmethod
    def _read_genomic_masked_regions(self) -> Dict[str, MaskedGenomicRegions]:
        pass
