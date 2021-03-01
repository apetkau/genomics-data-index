import abc
from pathlib import Path
from typing import List, Dict, Set

import pandas as pd


def verify_columns_match(expected_columns: Set[str], result_df: pd.DataFrame) -> None:
    actual_columns = set(result_df.columns.tolist())
    if not expected_columns.issubset(actual_columns):
        raise Exception('Results table does not contain expected set of columns. '
                        f'Expected {expected_columns}, actual {actual_columns}')


class QueryFeature:
    def __init__(self):
        pass


class QueryService(abc.ABC):

    def __init__(self):
        pass

    # Scenerio 1 (and 2 with proper implementation)
    def find_matches_genome_files(self, sample_reads: Dict[str, List[Path]],
                                  distance_threshold: float = None) -> pd.DataFrame:
        matches_df = self._find_matches_genome_files_internal(sample_reads, distance_threshold)
        matches_df.insert(loc=0, column='Type', value=self.get_data_type())
        verify_columns_match({'Type', 'Sample A', 'Sample B', 'Distance'}, matches_df)

        return matches_df

    @abc.abstractmethod
    def _find_matches_genome_files_internal(self, sample_reads: Dict[str, List[Path]],
                                            distance_threshold: float = None) -> pd.DataFrame:
        pass

    # TODO: Scenario 3 (by subtype) later

    # Scenario 4
    def find_matches(self, samples: List[str],
                     distance_threshold: float = None) -> pd.DataFrame:
        if samples is None or len(samples) == 0:
            raise Exception(f'Cannot find matches of empty samples list')

        matches_df = self._find_matches_internal(samples, distance_threshold)
        matches_df.insert(loc=0, column='Type', value=self.get_data_type())
        verify_columns_match({'Type', 'Sample A', 'Sample B', 'Distance'}, matches_df)

        return matches_df

    @abc.abstractmethod
    def _find_matches_internal(self, samples, distance_threshold: float):
        pass

    # Scenario 5
    def pairwise_distance(self, samples: List[str]) -> pd.DataFrame:
        if samples is None or len(samples) == 0:
            raise Exception(f'Cannot find pairwise distances of empty samples list')

        distance_df = self._pairwise_distance_internal(samples)
        distance_df.insert(loc=0, column='Type', value=self.get_data_type())
        verify_columns_match({'Type', 'Sample A', 'Sample B', 'Distance'}, distance_df)

        return distance_df

    @abc.abstractmethod
    def _pairwise_distance_internal(self, samples: List[str]) -> pd.DataFrame:
        pass

    # TODO: Scenario 6 later

    def find_by_features(self, features: List[QueryFeature]) -> pd.DataFrame:
        if features is None or len(features) == 0:
            raise Exception(f'Cannot find by empty features')

        matches_df = self._find_by_features_internal(features)
        matches_df.insert(loc=0, column='Type', value=self.get_data_type())
        verify_columns_match({'Type', 'Feature', 'Sample Name', 'Sample ID'}, matches_df)

        return matches_df

    def _find_by_features_internal(self, features: List[QueryFeature]) -> pd.DataFrame:
        pass

    # Scenario 8 in TreeService

    # Scenario 9
    def differences_between_genomes(self, sample1: str, sample2: str):
        differences_df = self._differences_between_genomes_internal(sample1, sample2)
        differences_df.insert(loc=0, column='Type', value=self.get_data_type())
        verify_columns_match({'Type', 'Sample1 unique', 'Sample2 unique'}, differences_df)

        return differences_df

    @abc.abstractmethod
    def _differences_between_genomes_internal(self, sample1: str, sample2: str):
        pass

    @abc.abstractmethod
    def get_data_type(self) -> str:
        pass
