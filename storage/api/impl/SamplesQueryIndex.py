from __future__ import annotations

from typing import Union

import pandas as pd

from storage.api.SamplesQuery import SamplesQuery
from storage.api.impl.QueriesCollection import QueriesCollection
from storage.connector.DataIndexConnection import DataIndexConnection
from storage.variant.SampleSet import SampleSet
from storage.variant.model.QueryFeature import QueryFeature
from storage.variant.model.QueryFeatureMLST import QueryFeatureMLST
from storage.variant.model.QueryFeatureMutation import QueryFeatureMutation
from storage.api.impl.TreeSamplesQuery import TreeSamplesQuery


class SamplesQueryIndex(SamplesQuery):

    HAS_KINDS = ['mutation', 'mlst']

    def __init__(self, connection: DataIndexConnection,
                 universe_set: SampleSet,
                 sample_set: SampleSet,
                 queries_collection: QueriesCollection = QueriesCollection.create_empty()):
        super().__init__()
        self._query_connection = connection
        self._universe_set = universe_set
        self._sample_set = sample_set
        self._queries_collection = queries_collection

    @property
    def universe_set(self) -> SampleSet:
        return self._universe_set

    @property
    def sample_set(self) -> SampleSet:
        return self._sample_set

    def intersect(self, sample_set: SampleSet, query_message: str = None) -> SamplesQuery:
        intersected_set = self._intersect_sample_set(sample_set)

        if query_message is None:
            query_message = f'intersect(samples={len(sample_set)}'

        queries_collection = self._queries_collection.append(query_message)
        return self._create_from(connection=self._query_connection, sample_set=intersected_set,
                                 queries_collection=queries_collection)

    def _intersect_sample_set(self, other: SampleSet) -> SampleSet:
        return self.sample_set.intersection(other)

    def build_tree(self, kind: str, scope: str, **kwargs):
        return TreeSamplesQuery.create(kind=kind, scope=scope, database_connection=self._query_connection,
                                       wrapped_query=self, **kwargs)

    def toframe(self, exclude_absent: bool = True) -> pd.DataFrame:
        sample_service = self._query_connection.sample_service
        return sample_service.create_dataframe_from_sample_set(self.sample_set,
                                                               universe_set=self._universe_set,
                                                               exclude_absent=exclude_absent,
                                                               queries_collection=self._queries_collection)

    def summary(self) -> pd.DataFrame:
        present = len(self)
        total = len(self._universe_set)
        unknown = pd.NA
        absent = total - present
        per_present = (present / total) * 100
        per_absent = (absent / total) * 100
        per_unknown = pd.NA

        return pd.DataFrame([{
            'Query': self._queries_collection.query_expression(),
            'Present': present,
            'Absent': absent,
            'Unknown': unknown,
            'Total': total,
            '% Present': per_present,
            '% Absent': per_absent,
            '% Unknown': per_unknown,
        }])

    def and_(self, other):
        if isinstance(other, SamplesQuery):
            intersect_set = self._intersect_sample_set(other.sample_set)
            queries_collection = self._queries_collection.append(str(other))
            return self._create_from(self._query_connection, intersect_set,
                                     queries_collection=queries_collection)
        else:
            raise Exception(f'Cannot perform an "and" on object {other}')

    def is_empty(self) -> bool:
        return self.sample_set.is_empty()

    def complement(self):
        complement_set = self.universe_set.minus(self.sample_set)
        query_collection = self._queries_collection.append('complement')
        return self._create_from(connection=self._query_connection,
                                 sample_set=complement_set,
                                 queries_collection=query_collection)

    @property
    def tree(self):
        raise Exception(f'No tree exists for {self.__class__}.'
                        f' Perhaps you should try to run build_tree() first to build a tree.')

    def __and__(self, other):
        return self.and_(other)

    def __len__(self):
        return len(self.sample_set)

    def __str__(self) -> str:
        percent_selected = (len(self) / len(self._universe_set)) * 100
        return f'<SamplesQueryIndex({len(self)}/{len(self._universe_set)} ({percent_selected:0.1f}) samples)>'

    def __repr__(self) -> str:
        return str(self)

    def has(self, feature: Union[QueryFeature, str], kind = None) -> SamplesQuery:
        if isinstance(feature, QueryFeature):
            query_feature = feature
        elif kind is None:
            raise Exception(f'feature=[{feature}] is not of type QueryFeature so must set "kind" parameter')
        elif kind == 'mutation':
            query_feature = QueryFeatureMutation(feature)
        elif kind == 'mlst':
            query_feature = QueryFeatureMLST(feature)
        else:
            raise Exception(f'kind={kind} is not recognized. Must be one of {self.HAS_KINDS}')

        found_set_dict = self._query_connection.sample_service.find_sample_sets_by_features([query_feature])

        if query_feature.id in found_set_dict:
            found_set = found_set_dict[query_feature.id]
            intersect_found = self._intersect_sample_set(found_set)
        else:
            intersect_found = SampleSet.create_empty()

        queries_collection = self._queries_collection.append(query_feature)
        return self._create_from(self._query_connection, intersect_found,
                                 queries_collection=queries_collection)

    def within(self, distance: float, sample_name: str, units: str) -> SamplesQuery:
        raise Exception(f'Cannot query within a distance without a tree.'
                        f' Perhaps you want to run build_tree() first to build a tree.')

    def is_type(self, sample_type) -> SamplesQuery:
        raise Exception('Not implemented')

    def _create_from(self, connection: DataIndexConnection, sample_set: SampleSet,
                     queries_collection: QueriesCollection) -> SamplesQuery:
        return SamplesQueryIndex(connection=self._query_connection,
                                 universe_set=self._universe_set,
                                 sample_set=sample_set,
                                 queries_collection=queries_collection)
