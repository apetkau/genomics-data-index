import re
from pathlib import Path

import pandas as pd

from genomics_data_index.storage.io.mlst.MLSTFeaturesReader import MLSTFeaturesReader
from genomics_data_index.storage.model import MLST_UNKNOWN_ALLELE


class MLSTTSeemannFeaturesReader(MLSTFeaturesReader):
    """
    A reader for results from the MLST software developed by Torsten Seemann (https://github.com/tseemann/mlst).
    Assumes output has been produced like:

    mlst --nopath *.fasta > mlst.tsv
    """

    def __init__(self, mlst_file: Path):
        super().__init__()

        self._mlst_file = mlst_file

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

    def _is_valid_allele(self, allele: str) -> bool:
        return allele != MLST_UNKNOWN_ALLELE and bool(re.match(r'^\d+$', allele))

    def _get_sample_from_filename(self, filename_series: pd.Series) -> pd.Series:
        file_sample_name_regex = r'^([^.]*)'
        return filename_series.str.extract(file_sample_name_regex, expand=True)

    def _extract_locus_alleles(self, df: pd.DataFrame) -> pd.DataFrame:
        locus_allele_list = list(set(df.columns) - {'File', 'Sample', 'Scheme', 'Sequence Type'})
        df = df.melt(id_vars=['File', 'Sample', 'Scheme', 'Sequence Type'], value_vars=locus_allele_list)
        df[['Locus', 'Allele']] = df['value'].str.extract(r'^([^\(]*)\(([^\)]*)\)', expand=True)
        return df
