import logging
from typing import Generator, Set

import pandas as pd

from genomics_data_index.storage.io.SampleData import SampleData
from genomics_data_index.storage.io.SampleDataPackage import SampleDataPackage
from genomics_data_index.storage.io.mlst.MLSTFeaturesReader import MLSTFeaturesReader
from genomics_data_index.storage.io.mlst.MLSTSampleData import MLSTSampleData

logger = logging.getLogger(__name__)


class MLSTSampleDataPackage(SampleDataPackage):

    def __init__(self, mlst_reader: MLSTFeaturesReader, index_unknown_missing: bool = True):
        super().__init__(index_unknown_missing=index_unknown_missing,
                         features_reader=mlst_reader)
        self._mlst_reader = mlst_reader

    def _minimal_expected_columns(self) -> Set[str]:
        return {'Sample', 'Scheme', 'Locus', 'Allele'}

    def _read_features_table(self) -> pd.DataFrame:
        return self._mlst_reader.get_features_table()

    def get_features_reader(self) -> MLSTFeaturesReader:
        if not self.index_unknown_missing():
            logger.warning(f'index_unknown_missing=False for MLSTSampleDataPackage, but right now '
                           f'missing data for MLST results is always indexed so this will be ignored.')
        return self._mlst_reader

    def sample_names(self) -> Set[str]:
        return set(self._mlst_reader.samples_list())

    def iter_sample_data(self) -> Generator[SampleData, None, None]:
        for sample_name, scheme, data in self._mlst_reader.iter_sample_data():
            yield MLSTSampleData(sample_name=sample_name, scheme=scheme, mlst_allele_data=data)
