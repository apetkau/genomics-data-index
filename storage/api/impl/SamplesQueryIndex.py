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


class SamplesQueryIndex(SamplesQuery):

    HAS_KINDS = ['mutation', 'mlst']

    def __init__(self, connection: DataIndexConnection,
                 sample_set: SampleSet = SamplesQuery.ALL_SAMPLES,
                 queries_collection: QueriesCollection = QueriesCollection.create_empty()):
        super().__init__()
        self._query_connection = connection
        self._sample_set = sample_set
        self._queries_collection = queries_collection

    @property
    def sample_set(self) -> SampleSet:
        return self._sample_set

    def _intersect_sample_set(self, other: SampleSet) -> SampleSet:
        return self.sample_set.intersection(other)

    def toframe(self) -> pd.DataFrame:
        sample_service = self._query_connection.sample_service
        return sample_service.create_dataframe_from_sample_set(self.sample_set,
                                                               self._queries_collection)

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

    def __and__(self, other):
        return self.and_(other)

    def __len__(self):
        return len(self.sample_set)

    def __str__(self) -> str:
        return f'<SamplesQueryIndex(len={len(self)})>'

    def __repr__(self) -> str:
        return f'<SamplesQueryIndex(len={len(self)})>'

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

    def within(self, distance: float, sample_name: str, kind: str) -> SamplesQuery:
        raise Exception('Not implemented')

    def within_kmer(self, distance: float, sample_name: str) -> SamplesQuery:
        return self.within(distance, sample_name, 'kmer')

    def is_type(self, sample_type) -> SamplesQuery:
        raise Exception('Not implemented')

    def _create_from(self, connection: DataIndexConnection, sample_set: SampleSet,
                     queries_collection: QueriesCollection) -> SamplesQuery:
        return SamplesQueryIndex(connection=self._query_connection,
                                 sample_set=sample_set,
                                 queries_collection=queries_collection)
