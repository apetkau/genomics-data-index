import abc
from typing import Dict

import pandas as pd

from storage.variant.io.FeaturesReader import FeaturesReader
from storage.variant.MaskedGenomicRegions import MaskedGenomicRegions


class NucleotideFeaturesReader(FeaturesReader):

    def __init(self):
        super().__init__()

    def _check_features_table_columns(self, features_df: pd.DataFrame) -> None:
        expected_columns = {'SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE'}
        actual_columns = set(features_df.columns.tolist())
        if not expected_columns.issubset(actual_columns):
            raise Exception('Variants table does not contain expected set of columns. '
                            f'Expected {expected_columns}, actual {actual_columns}')

    def get_genomic_masked_regions(self) -> Dict[str, MaskedGenomicRegions]:
        return self._read_genomic_masked_regions()

    @abc.abstractmethod
    def _read_genomic_masked_regions(self) -> Dict[str, MaskedGenomicRegions]:
        pass
