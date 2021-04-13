from pathlib import Path
from typing import Dict, List

import pandas as pd

from storage.variant.io.MLSTFeaturesReader import MLSTFeaturesReader


class BasicMLSTFeaturesReader(MLSTFeaturesReader):

    def __init__(self, mlst_file: Path):
        super().__init__()

        self._mlst_file = mlst_file
        self._features_table = None

    def sample_feature_files(self) -> Dict[str, Path]:
        raise Exception('Not implemented')

    def samples_list(self) -> List[str]:
        mlst_df = self.get_features_table()
        return list(set(mlst_df['Sample'].tolist()))

    def get_features_table(self) -> pd.DataFrame:
        if self._features_table is None:
            self._features_table = super().get_features_table()
        return self._features_table

    def _read_features_table(self) -> pd.DataFrame:
        df = pd.read_csv(self._mlst_file, sep='\t', header=None)
        df = df.rename(columns={
            0: 'File',
            1: 'Scheme',
            2: 'Sequence Type',
        })

        df['Sample'] = self._get_sample_from_filename(df['File'])
        df = self._extract_locus_alleles(df)

        df = df[['File', 'Sample', 'Scheme', 'Locus', 'Allele', 'Sequence Type']].sort_values(
            by=['Sample', 'Scheme', 'Locus']).reset_index().drop(columns='index')

        return df

    def _get_sample_from_filename(self, filename_series: pd.Series) -> pd.Series:
        file_sample_name_regex = r'^([^.]*)'
        return filename_series.str.extract(file_sample_name_regex, expand=True)

    def _extract_locus_alleles(self, df: pd.DataFrame) -> pd.DataFrame:
        locus_allele_list = list(set(df.columns) - {'File', 'Sample', 'Scheme', 'Sequence Type'})
        df = df.melt(id_vars=['File', 'Sample', 'Scheme', 'Sequence Type'], value_vars=locus_allele_list)
        df[['Locus', 'Allele']] = df['value'].str.extract(r'^([^\(]*)\(([^\)]*)\)', expand=True)
        return df
