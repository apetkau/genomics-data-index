from __future__ import annotations

import abc
from typing import Union, List, Set, Tuple

import numpy as np
import pandas as pd
from ete3 import Tree

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
             sample_names_column: str = None, default_isa_kind: str = 'names',
             default_isa_column: str = None) -> SamplesQuery:
        """
        Joins the passed dataframe onto the current query using the passed column name in the dataframe.
        The column can either contain sample IDs (in sample_ids_column) or sample names (sample_names_column).
        This will modify the universe set to the subset of samples found within the passed data frame.
        :param data_frame: The data frame to join on.
        :param sample_ids_column: The column name in the data frame containing internal sample IDs used for joining.
        :param sample_names_column: The column name in the data frame containing sample names to join on.
                                    Internally these will be mapped to Sample IDs.
        :param default_isa_kind: The default 'kind' parameter to isa(data, kind=kind, ...) queries. This is used
                                 to simplify queries after joining with a dataframe in cases where you want
                                 isa() to select by a particular column in the dataframe.
        :param default_isa_column: If default_isa_kind is set to 'dataframe' sets the default column to use for
                                   isa() queries.
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
    def or_(self, other: SamplesQuery) -> SamplesQuery:
        pass

    @abc.abstractmethod
    def build_tree(self, kind: str, **kwargs) -> SamplesQuery:
        pass

    @abc.abstractmethod
    def join_tree(self, tree: Tree, kind='mutation', **kwargs) -> SamplesQuery:
        pass

    @abc.abstractmethod
    def reset_universe(self) -> SamplesQuery:
        """
        Resets the 'universe' set to be the set of currently selected samples.
        :return: A SamplesQuery with the universe reset to whatever is selected.
        """
        pass

    @abc.abstractmethod
    def hasa(self, property: Union[QueryFeature, str, pd.Series], kind='mutation') -> SamplesQuery:
        pass

    def has(self, property: Union[QueryFeature, str, pd.Series], kind='mutation') -> SamplesQuery:
        """
        Queries for samples that have a particular property. Synonym for hasa().
        """
        return self.hasa(property=property, kind=kind)

    @abc.abstractmethod
    def complement(self):
        pass

    def within(self, data: Union[str, List[str], SamplesQuery, SampleSet], **kwargs) -> SamplesQuery:
        """
        Queries for samples within a particular distance. This is identical to calling
        isin(data, kind='distance', ...). Used to help make code easier to read.

        :param data: The data to use for selecting samples by.
        :param **kwargs: Other arguments for the the internal implementations.
        :return: A SamplesQuery with the matched samples.
        """
        return self.isin(data=data, kind='distance', **kwargs)

    @abc.abstractmethod
    def _within_distance(self, data: Union[str, List[str], SamplesQuery, SampleSet], distance: float, units: str,
                         **kwargs) -> SamplesQuery:
        pass

    @abc.abstractmethod
    def _distance_units(self) -> List[str]:
        pass

    @abc.abstractmethod
    def isa(self, data: Union[str, List[str], SamplesQuery, SampleSet], kind: str = None,
            **kwargs) -> SamplesQuery:
        """
        Queries for samples which are a particular type/belong to a particular category.
        Read as "subset samples which are a (isa) particular type defined by 'data'".
        The default implementation will select samples by sample name but this is useful when used
        with an attached dataframe (at which point it selects based on matches to a column,
        see documentation for DataFrameSamplesQuery).

        :param data: The data to match.
        :param kind: The particular kind of data passed.
        :**kwargs: Arguments for the underlying implementation.
        :return: A SamplesQuery with the matched samples.
        """
        pass

    def isan(self, data: Union[str, List[str], SamplesQuery, SampleSet], kind: str = None,
             **kwargs) -> SamplesQuery:
        """
        Synonym for isa()
        """
        return self.isa(data=data, kind=kind, **kwargs)

    @abc.abstractmethod
    def isin(self, data: Union[str, List[str], pd.Series, SamplesQuery, SampleSet], kind: str = None,
             **kwargs) -> SamplesQuery:
        pass

    @abc.abstractmethod
    def _get_sample_names_query_infix_from_data(self, data: Union[str, List[str], pd.Series, SamplesQuery, SampleSet]
                                            ) -> Tuple[Set[str], str]:
        pass

    @abc.abstractmethod
    def _can_handle_isin_kind(self, kind: str) -> bool:
        pass

    @abc.abstractmethod
    def _can_handle_isa_kind(self, kind: str) -> bool:
        pass

    @abc.abstractmethod
    def _isa_kinds(self) -> List[str]:
        pass

    @abc.abstractmethod
    def _isin_kinds(self) -> List[str]:
        pass

    @abc.abstractmethod
    def to_distances(self, kind: str, **kwargs) -> Tuple[np.ndarray, List[str]]:
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
    def has_tree(self) -> bool:
        pass

    @abc.abstractmethod
    def tolist(self, names=True) -> List[Union[int, str]]:
        pass

    @abc.abstractmethod
    def __and__(self, other):
        pass

    @abc.abstractmethod
    def __or__(self, other):
        pass

    @abc.abstractmethod
    def __len__(self):
        pass
