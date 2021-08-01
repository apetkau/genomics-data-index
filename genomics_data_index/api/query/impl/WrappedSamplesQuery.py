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

    def __init__(self, connection: DataIndexConnection, wrapped_query: SamplesQuery):
        """
        Builds a new WrappedSamplesQuery from the given information. In most normal operations WrappedSamplesQuery objects
        are not created directly but are instead created from an :py:class:`genomics_data_index.api.GenomicsDataIndex`
        object or from operations applied to a SamplesQuery.

        :param connection: A connection to a database containing samples.
        :param wrapped_query: The SamplesQuery to wrap around/decorate.
        :return: A new WrappedSamplesQuery object.
        """
        super().__init__()
        self._query_connection = connection
        self._wrapped_query = wrapped_query

    @property
    def universe_set(self) -> SampleSet:
        return self._wrapped_query.universe_set

    @property
    def sample_set(self) -> SampleSet:
        return self._wrapped_query.sample_set

    @property
    def unknown_set(self) -> SampleSet:
        return self._wrapped_query.unknown_set

    @property
    def absent_set(self) -> SampleSet:
        return self._wrapped_query.absent_set

    def select_absent(self) -> SamplesQuery:
        return self._wrap_create(self._wrapped_query.select_absent())

    def select_unknown(self) -> SamplesQuery:
        return self._wrap_create(self._wrapped_query.select_unknown())

    def select_present(self) -> SamplesQuery:
        return self._wrap_create(self._wrapped_query.select_present())

    def reset_universe(self, include_unknown: bool = True) -> SamplesQuery:
        return self._wrap_create(self._wrapped_query.reset_universe())

    def set_universe(self, universe_set: SampleSet, query_message: str = None) -> SamplesQuery:
        return self._wrap_create(self._wrapped_query.set_universe(universe_set, query_message=query_message))

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
    def _wrap_create(self, wrapped_query: SamplesQuery) -> WrappedSamplesQuery:
        pass

    def toframe(self, include_present: bool = True, include_unknown: bool = False,
                include_absent: bool = False) -> pd.DataFrame:
        return self._wrapped_query.toframe(include_present=include_present, include_unknown=include_unknown,
                                           include_absent=include_absent)

    def summary(self) -> pd.DataFrame:
        return self._wrapped_query.summary()

    def features_summary(self, kind: str = 'mutations', selection: str = 'all',
                         include_present_features: bool = True, include_unknown_features: bool = False,
                         **kwargs) -> pd.DataFrame:
        return self._wrapped_query.features_summary(kind=kind, selection=selection,
                                                    include_present_features=include_present_features,
                                                    include_unknown_features=include_unknown_features,
                                                    **kwargs)

    def tofeaturesset(self, kind: str = 'mutations', selection: str = 'all',
                      include_present_features: bool = True, include_unknown_features: bool = False) -> Set[str]:
        return self._wrapped_query.tofeaturesset(kind=kind, selection=selection,
                                                 include_present_features=include_present_features,
                                                 include_unknown_features=include_unknown_features)

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

    def tolist(self, names: bool = True, include_present: bool = True,
               include_unknown: bool = False, include_absent: bool = False) -> Union[List[str], List[int]]:
        return self._wrapped_query.tolist(names=names, include_present=include_present,
                                          include_unknown=include_unknown, include_absent=include_absent)

    def toset(self, names: bool = True, include_present: bool = True,
              include_unknown: bool = False, include_absent: bool = False) -> Union[Set[str], Set[int]]:
        return self._wrapped_query.toset(names=names, include_present=include_present,
                                         include_unknown=include_unknown, include_absent=include_absent)

    def _get_sample_name_ids(self, include_unknowns: bool = True) -> Dict[str, int]:
        sample_service = self._query_connection.sample_service
        name_ids = {s.name: s.id for s in sample_service.find_samples_by_ids(self._wrapped_query.sample_set)}

        if include_unknowns:
            name_ids_unknown = {s.name: s.id for s in
                                sample_service.find_samples_by_ids(self._wrapped_query.unknown_set)}
            name_ids.update(name_ids_unknown)

        return name_ids

    def is_empty(self, include_unknown=False) -> bool:
        return self._wrapped_query.is_empty(include_unknown=include_unknown)

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

    def __len__(self):
        return len(self._wrapped_query)

    def __str__(self) -> str:
        universe_length = len(self.universe_set)
        unknown_length = len(self.unknown_set)
        selected_length = len(self)

        if universe_length > 0:
            percent_selected = f'{(selected_length / universe_length) * 100:0.0f}'
            percent_unknown = f'{(unknown_length / universe_length) * 100:0.0f}'
        else:
            percent_selected = 'NA'
            percent_unknown = 'NA'

        return f'<{self.__class__.__name__}[selected={percent_selected}% ' \
               f'({selected_length}/{universe_length}) samples, unknown={percent_unknown}% ' \
               f'({unknown_length}/{universe_length}) samples]>'

    def __repr__(self) -> str:
        return str(self)
