from pathlib import Path

import pandas as pd

from genomics_data_index.storage.io.SampleData import SampleData


class MLSTSampleData(SampleData):

    def __init__(self, sample_name: str, scheme: str, mlst_allele_data: pd.Series):
        super().__init__(sample_name=sample_name)
        self._mlst_allele_data = mlst_allele_data
        self._scheme = scheme

    def get_scheme(self) -> str:
        return self._scheme

    def get_alleles(self):
        return self._mlst_allele_data

    def is_preprocessed(self) -> bool:
        return True

    def _do_preprocess_and_persist(self, output_dir: Path) -> SampleData:
        return self

    def _do_persist(self, output_dir: Path) -> SampleData:
        return self
