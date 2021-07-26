import abc
from typing import Generator, Set, Dict
import pandas as pd

from genomics_data_index.storage.io.SampleData import SampleData
from genomics_data_index.storage.io.FeaturesReader import FeaturesReader


class SampleDataPackage(abc.ABC):

    def __init__(self, index_unknown_missing: bool,
                 features_reader: FeaturesReader):
        self._index_unknown_missing = index_unknown_missing
        self._features_reader = features_reader

    @abc.abstractmethod
    def sample_names(self) -> Set[str]:
        pass

    def index_unknown_missing(self) -> bool:
        return self._index_unknown_missing

    def _check_features_table_columns(self, features_df: pd.DataFrame) -> None:
        expected_columns = self._minimal_expected_columns()
        actual_columns = set(features_df.columns.tolist())
        if not expected_columns.issubset(actual_columns):
            raise Exception('Variants table does not contain expected set of columns. '
                            f'Expected {expected_columns}, actual {actual_columns}')

    @abc.abstractmethod
    def _minimal_expected_columns(self) -> Set[str]:
        pass

    @abc.abstractmethod
    def _read_features_table(self) -> pd.DataFrame:
        pass

    def get_features_reader(self):
        return self._features_reader

    def get_features_table(self) -> pd.DataFrame:
        features_df = self._read_features_table()
        self._check_features_table_columns(features_df)
        return features_df

    def process_all_data(self) -> Dict[str, SampleData]:
        processed_data = {}
        for sample_data in self.iter_sample_data():
            processed_data[sample_data.sample_name] = sample_data

        return processed_data

    @abc.abstractmethod
    def iter_sample_data(self) -> Generator[SampleData, None, None]:
        pass
