from typing import Dict, List, Set
import abc

import pandas as pd
from pathlib import Path


class FeaturesReader(abc.ABC):

    def __init(self):
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
    def sample_feature_files(self) -> Dict[str, Path]:
        """
        Gets a dictionary of sample names to feature files to be read by this reader.
        :return: A dictionary of sample names to feature files ('name' => 'file')
        """
        pass

    @abc.abstractmethod
    def samples_list(self) -> List[str]:
        """
        Gets a list of sample names that will be read by this reader.
        :return: A list of sample names that will be read by this reader.
        """
        pass

    @abc.abstractmethod
    def _read_features_table(self) -> pd.DataFrame:
        pass
