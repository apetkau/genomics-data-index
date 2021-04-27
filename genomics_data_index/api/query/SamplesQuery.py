from __future__ import annotations

import abc
from typing import Union, List, Set

import pandas as pd

from genomics_data_index.storage.SampleSet import SampleSet
from genomics_data_index.storage.model.QueryFeature import QueryFeature


class SamplesQuery(abc.ABC):
    ALL_SAMPLES = SampleSet.create_all()

    def __init__(self):
        pass

    @property
    @abc.abstractmethod
    def universe_set(self) -> SampleSet:
        pass

    @property
    @abc.abstractmethod
    def sample_set(self) -> SampleSet:
        pass

    @abc.abstractmethod
    def join(self, data_frame: pd.DataFrame, sample_ids_column: str = None,
             sample_names_column: str = None) -> SamplesQuery:
        """
        Joins the passed dataframe onto the current query using the passed column name in the dataframe.
        The column can either contain sample IDs (in sample_ids_column) or sample names (sample_names_column).
        This will modify the universe set to the subset of samples found within the passed data frame.
        :param data_frame: The data frame to join on.
        :param sample_ids_column: The column name in the data frame containing internal sample IDs used for joining.
        :param sample_names_column: The column name in the data frame containing sample names to join on.
                                    Internally these will be mapped to Sample IDs.
        :return: A SamplesQuery representing the query joined with the data frame.
        """
        pass

    @abc.abstractmethod
    def toframe(self, exclude_absent: bool = True) -> pd.DataFrame:
        pass

    @abc.abstractmethod
    def summary(self) -> pd.DataFrame:
        pass

    @abc.abstractmethod
    def summary_features(self, kind: str = 'mutations', **kwargs) -> pd.DataFrame:
        pass

    @abc.abstractmethod
    def tofeaturesset(self, kind: str = 'mutations', selection: str = 'all',
                      ncores: int = 1) -> Set[str]:
        pass

    @abc.abstractmethod
    def intersect(self, sample_set: SampleSet, query_message: str = None) -> SamplesQuery:
        pass

    @abc.abstractmethod
    def and_(self, other: SamplesQuery) -> SamplesQuery:
        pass

    @abc.abstractmethod
    def build_tree(self, kind: str, scope: str, **kwargs) -> SamplesQuery:
        pass

    @abc.abstractmethod
    def hasa(self, property: Union[QueryFeature, str, pd.Series], kind=None) -> SamplesQuery:
        pass

    def has(self, property: Union[QueryFeature, str, pd.Series], kind=None) -> SamplesQuery:
        return self.hasa(property=property, kind=kind)

    @abc.abstractmethod
    def complement(self):
        pass
    #
    # @abc.abstractmethod
    # def within(self, sample_names: Union[str, List[str]], kind: str = 'distance', **kwargs) -> SamplesQuery:
    #     return self.isin(sample_names=sample_names, kind='distance')

    @abc.abstractmethod
    def isin(self, sample_names: Union[str, List[str]], kind: str = 'distance', **kwargs) -> SamplesQuery:
        pass

    @abc.abstractmethod
    def _get_has_kinds(self) -> List[str]:
        pass

    @abc.abstractmethod
    def is_empty(self):
        pass

    @abc.abstractmethod
    def query_expression(self) -> str:
        pass

    @property
    @abc.abstractmethod
    def tree(self):
        pass

    @abc.abstractmethod
    def tolist(self, names=True) -> List[Union[int, str]]:
        pass

    @abc.abstractmethod
    def __and__(self, other):
        pass

    @abc.abstractmethod
    def __len__(self):
        pass
