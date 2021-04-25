import re
from pathlib import Path

import pandas as pd

from genomics_data_index.storage.io.mlst.MLSTFeaturesReader import MLSTFeaturesReader
from genomics_data_index.storage.model import MLST_UNKNOWN_ALLELE


class MLSTSistrReader(MLSTFeaturesReader):
    SCHEME_NAME = 'sistr_330'

    def __init__(self, mlst_file: Path):
        super().__init__()

        self._mlst_file = mlst_file

    def _read_features_table(self) -> pd.DataFrame:
        df = pd.read_csv(self._mlst_file, sep=',', dtype=str)
        df = df.rename(columns={
            'Unnamed: 0': 'File',
        })

        df['Sample'] = self._get_sample_from_filename(df['File'])
        df = self._extract_locus_alleles(df)
        df['Scheme'] = self.SCHEME_NAME

        df = df[['File', 'Sample', 'Scheme', 'Locus', 'Allele']].sort_values(
            by=['Sample', 'Scheme', 'Locus']).reset_index().drop(columns='index')

        return df

    def _is_valid_allele(self, allele: str) -> bool:
        return allele != MLST_UNKNOWN_ALLELE and bool(re.match(r'^\d+$', allele))

    def _get_sample_from_filename(self, filename_series: pd.Series) -> pd.Series:
        file_sample_name_regex = r'^([^.]*)'
        return filename_series.str.extract(file_sample_name_regex, expand=True)

    def _extract_locus_alleles(self, df: pd.DataFrame) -> pd.DataFrame:
        locus_allele_list = list(set(df.columns) - {'File', 'Sample'})
        df = df.melt(id_vars=['File', 'Sample'], value_vars=locus_allele_list,
                     var_name='Locus', value_name='Allele')
        df['Allele'] = df['Allele'].astype(str)
        return df
