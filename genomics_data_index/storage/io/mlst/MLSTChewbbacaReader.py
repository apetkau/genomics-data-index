import re
from pathlib import Path

import pandas as pd

from genomics_data_index.storage.io.mlst.MLSTFeaturesReader import MLSTFeaturesReader
from genomics_data_index.storage.model import MLST_UNKNOWN_ALLELE


class MLSTChewbbacaReader(MLSTFeaturesReader):

    def __init__(self, mlst_file: Path, scheme: str):
        super().__init__()

        self._mlst_file = mlst_file

        if scheme is None:
            raise Exception('scheme cannot be None')

        self._scheme = scheme

    def _read_features_table(self) -> pd.DataFrame:
        df = pd.read_csv(self._mlst_file, sep='\t', dtype=str)
        df = df.rename(columns={
            'FILE': 'File',
        })

        df['Sample'] = self._get_sample_from_filename(df['File'])
        df = self._extract_locus_alleles(df)
        df['Scheme'] = self._scheme

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
