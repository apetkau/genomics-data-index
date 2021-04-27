from __future__ import annotations

import abc
from typing import Set, Union, Dict, List

import pandas as pd

from genomics_data_index.api.query.SamplesQuery import SamplesQuery
from genomics_data_index.configuration.connector import DataIndexConnection
from genomics_data_index.storage.SampleSet import SampleSet
from genomics_data_index.storage.model.QueryFeature import QueryFeature


class WrappedSamplesQuery(SamplesQuery, abc.ABC):

    def __init__(self, connection: DataIndexConnection, wrapped_query: SamplesQuery, universe_set: SampleSet = None):
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

    def intersect(self, sample_set: SampleSet, query_message: str = None) -> SamplesQuery:
        intersected_query = self._wrapped_query.intersect(sample_set=sample_set, query_message=query_message)
        return self._wrap_create(intersected_query)

    def join(self, data_frame: pd.DataFrame, sample_ids_column: str = None,
             sample_names_column: str = None) -> SamplesQuery:
        return self._wrap_create(self._wrapped_query.join(data_frame=data_frame, sample_ids_column=sample_ids_column,
                                                          sample_names_column=sample_names_column))

    @abc.abstractmethod
    def _wrap_create(self, wrapped_query: SamplesQuery) -> WrappedSamplesQuery:
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

    def hasa(self, property: Union[QueryFeature, str, pd.Series], kind=None) -> SamplesQuery:
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

    def is_type(self, sample_type) -> SamplesQuery:
        return self._wrap_create(self._wrapped_query.is_type(sample_type))

    def is_empty(self):
        return self._wrapped_query.is_empty()

    @property
    def tree(self):
        raise Exception(f'No tree exists for {self.__class__}.'
                        f' Perhaps you should try to run build_tree() first to build a tree.')

    def isin(self, sample_names: Union[str, List[str]], kind: str = 'distance', **kwargs) -> SamplesQuery:
        raise Exception(f'Cannot query within a distance without a tree.'
                        f' Perhaps you want to run build_tree() first to build a tree.')

    def __and__(self, other):
        return self.and_(other)

    def __len__(self):
        return len(self._wrapped_query)

    def __str__(self) -> str:
        universe_length = len(self.universe_set)
        percent_selected = (len(self) / universe_length) * 100
        return f'<{self.__class__.__name__}[{percent_selected:0.0f}% ({len(self)}/{universe_length}) samples]>'

    def __repr__(self) -> str:
        return str(self)
