import abc
from typing import List, Set, Optional

import pandas as pd

from genomics_data_index.storage.io.SampleData import SampleData


class FeaturesReader(abc.ABC):

    def __init__(self):
        pass

    def _check_features_table_columns(self, features_df: pd.DataFrame) -> None:
        expected_columns = self._minimal_expected_columns()
        actual_columns = set(features_df.columns.tolist())
        if not expected_columns.issubset(actual_columns):
            raise Exception('Variants table does not contain expected set of columns. '
                            f'Expected {expected_columns}, actual {actual_columns}')

    @abc.abstractmethod
    def _minimal_expected_columns(self) -> Set[str]:
        pass

    def get_features_table(self) -> pd.DataFrame:
        features_df = self._read_features_table()
        self._check_features_table_columns(features_df)
        return features_df

    @abc.abstractmethod
    def get_sample_files(self, sample_name: str) -> Optional[SampleData]:
        pass

    @abc.abstractmethod
    def get_or_create_feature_file(self, sample_name: str):
        """
        Gets a file of the features associated with the sample (or creates such a file if it doesn't exist).
        :param sample_name: The sample to get features for.
        :return: A Path to the file containing all the features for the sample.
        """
        pass

    @abc.abstractmethod
    def samples_list(self) -> List[str]:
        """
        Gets a list of sample names that will be read by this reader.
        :return: A list of sample names that will be read by this reader.
        """
        pass

    def samples_set(self) -> Set[str]:
        """
        Gets a set of sample names that will be read by this reader.
        :return: A set of sample names that will be read by this reader.
        """
        return set(self.samples_list())

    @abc.abstractmethod
    def _read_features_table(self) -> pd.DataFrame:
        pass
