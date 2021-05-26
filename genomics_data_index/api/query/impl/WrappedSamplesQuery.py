from __future__ import annotations

import abc
from typing import Set, Union, Dict, List, Tuple

import numpy as np
import pandas as pd

from genomics_data_index.api.query.SamplesQuery import SamplesQuery
from genomics_data_index.configuration.connector import DataIndexConnection
from genomics_data_index.storage.SampleSet import SampleSet
from genomics_data_index.storage.model.QueryFeature import QueryFeature


class WrappedSamplesQuery(SamplesQuery, abc.ABC):
    """
    A WrappedSamplesQuery wraps around/decorates a SamplesQuery object with additional information used to enhance the
    queries. For example by including a tree or a DataFrame.
    """

    def __init__(self, connection: DataIndexConnection, wrapped_query: SamplesQuery, universe_set: SampleSet = None):
        """
        Builds a new WrappedSamplesQuery from the given information. In most normal operations WrappedSamplesQuery objects
        are not created directly but are instead created from an :py:class:`genomics_data_index.api.GenomicsDataIndex`
        object or from operations applied to a SamplesQuery.

        :param connection: A connection to a database containing samples.
        :param wrapped_query: The SamplesQuery to wrap around/decorate.
        :param universe_set: The :py:class:`genomics_data_index.storage.SampleSet` representing a set of samples defining
                       the universe (used for e.g., complement() operations).
        :return: A new WrappedSamplesQuery object.
        """
        super().__init__()
        self._query_connection = connection
        self._wrapped_query = wrapped_query
        self._universe_set = universe_set

    @property
    def universe_set(self) -> SampleSet:
        if self._universe_set is None:
            return self._wrapped_query.universe_set
        else:
            return self._universe_set

    @property
    def sample_set(self) -> SampleSet:
        return self._wrapped_query.sample_set

    def reset_universe(self) -> SamplesQuery:
        return self._wrap_create(self._wrapped_query.reset_universe())

    def intersect(self, sample_set: SampleSet, query_message: str = None) -> SamplesQuery:
        intersected_query = self._wrapped_query.intersect(sample_set=sample_set, query_message=query_message)
        return self._wrap_create(intersected_query)

    def join(self, data_frame: pd.DataFrame, sample_ids_column: str = None,
             sample_names_column: str = None, default_isa_kind: str = 'names',
             default_isa_column: str = None) -> SamplesQuery:
        return self._wrap_create(self._wrapped_query.join(data_frame=data_frame, sample_ids_column=sample_ids_column,
                                                          sample_names_column=sample_names_column,
                                                          default_isa_kind=default_isa_kind,
                                                          default_isa_column=default_isa_column))

    @abc.abstractmethod
    def _wrap_create(self, wrapped_query: SamplesQuery, universe_set: SampleSet = None) -> WrappedSamplesQuery:
        pass

    def toframe(self, exclude_absent: bool = True) -> pd.DataFrame:
        return self._wrapped_query.toframe(exclude_absent=exclude_absent)

    def summary(self) -> pd.DataFrame:
        return self._wrapped_query.summary()

    def summary_features(self, kind: str = 'mutations', **kwargs) -> pd.DataFrame:
        return self._wrapped_query.summary_features(kind=kind, **kwargs)

    def tofeaturesset(self, kind: str = 'mutations', selection: str = 'all', ncores: int = 1) -> Set[str]:
        return self._wrapped_query.tofeaturesset(kind=kind, selection=selection, ncores=ncores)

    def and_(self, other: SamplesQuery) -> SamplesQuery:
        return self._wrap_create(self._wrapped_query.and_(other))

    def or_(self, other: SamplesQuery) -> SamplesQuery:
        return self._wrap_create(self._wrapped_query.or_(other))

    def hasa(self, property: Union[QueryFeature, str, pd.Series], kind='mutation') -> SamplesQuery:
        return self._wrap_create(self._wrapped_query.hasa(property=property, kind=kind))

    def _get_has_kinds(self) -> List[str]:
        return self._wrapped_query._get_has_kinds()

    def complement(self):
        return self._wrap_create(self._wrapped_query.complement())

    def query_expression(self) -> str:
        return self._wrapped_query.query_expression()

    def tolist(self, names=True):
        return self._wrapped_query.tolist(names=names)

    def _get_sample_name_ids(self) -> Dict[str, int]:
        sample_service = self._query_connection.sample_service
        return {s.name: s.id for s in sample_service.find_samples_by_ids(self._wrapped_query.sample_set)}

    def is_empty(self):
        return self._wrapped_query.is_empty()

    def has_tree(self) -> bool:
        return self._wrapped_query.has_tree()

    @property
    def tree(self):
        raise Exception(f'No tree exists for {self.__class__}.'
                        f' Perhaps you should try to run build_tree() first to build a tree.')

    def _isin_internal(self, data: Union[str, List[str], pd.Series], kind: str, **kwargs) -> SamplesQuery:
        return self._wrap_create(self._wrapped_query.isin(data=data, kind=kind, **kwargs))

    def _isa_internal(self, data: Union[str, List[str]], kind: str, **kwargs) -> SamplesQuery:
        return self._wrap_create(self._wrapped_query.isa(data=data, kind=kind, **kwargs))

    def _isin_kinds(self) -> List[str]:
        return list(set(self._wrapped_query._isin_kinds()))

    def _isa_kinds(self) -> List[str]:
        return list(set(self._wrapped_query._isa_kinds()))

    def isin(self, data: Union[str, List[str], SamplesQuery, SampleSet], kind: str = 'samples',
             **kwargs) -> SamplesQuery:
        if self._can_handle_isin_kind(kind):
            return self._isin_internal(data=data, kind=kind, **kwargs)
        else:
            return self._wrap_create(self._wrapped_query.isin(data=data, kind=kind, **kwargs))

    def _get_sample_names_query_infix_from_data(self, data: Union[str, List[str], pd.Series, SamplesQuery, SampleSet]
                                                ) -> Tuple[Set[str], str]:
        return self._wrapped_query._get_sample_names_query_infix_from_data(data)

    def _distance_units(self) -> List[str]:
        return self._wrapped_query._distance_units()

    def _within_distance(self, data: Union[str, List[str], SamplesQuery, SampleSet], distance: float,
                         units: str, **kwargs) -> SamplesQuery:
        return self._wrap_create(self._wrapped_query._within_distance(data=data,
                                                                      distance=distance,
                                                                      units=units,
                                                                      **kwargs))

    def _can_handle_isa_kind(self, kind: str) -> bool:
        return False

    def isa(self, data: Union[str, List[str]], kind: str = 'sample', **kwargs) -> SamplesQuery:
        if self._can_handle_isa_kind(kind):
            return self._isa_internal(data=data, kind=kind, **kwargs)
        else:
            return self._wrap_create(self._wrapped_query.isa(data=data, kind=kind, **kwargs))

    def to_distances(self, kind: str = 'kmer', **kwargs) -> Tuple[np.ndarray, List[str]]:
        return self._wrapped_query.to_distances(kind=kind, **kwargs)

    def __and__(self, other):
        return self.and_(other)

    def __or__(self, other):
        return self.or_(other)

    def __len__(self):
        return len(self._wrapped_query)

    def __str__(self) -> str:
        universe_length = len(self.universe_set)
        percent_selected = (len(self) / universe_length) * 100
        return f'<{self.__class__.__name__}[{percent_selected:0.0f}% ({len(self)}/{universe_length}) samples]>'

    def __repr__(self) -> str:
        return str(self)
