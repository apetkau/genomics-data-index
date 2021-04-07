from pathlib import Path
from typing import Dict, List

import pandas as pd

from storage.variant.io.MLSTFeaturesReader import MLSTFeaturesReader


class BasicMLSTFeaturesReader(MLSTFeaturesReader):

    def __init__(self, mlst_file: Path):
        super().__init__()

        self._mlst_file = mlst_file

    def sample_feature_files(self) -> Dict[str, Path]:
        raise Exception('Not implemented')

    def samples_list(self) -> List[str]:
        pass

    def _read_features_table(self) -> pd.DataFrame:
        df = pd.read_csv(self._mlst_file, sep='\t', header=None)
        df = df.rename(columns={
            0: 'File',
            1: 'Scheme',
            2: 'Sequence Type',
        })
        return df

    def add_sample_column(self, df: pd.DataFrame) -> pd.DataFrame:
        df['']
