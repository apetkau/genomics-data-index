from __future__ import annotations

import pandas as pd

from storage.variant.model.QueryFeature import QueryFeature
from storage.variant.model.QueryFeatureMutation import QueryFeatureMutation
from storage.variant.model.QueryFeatureMLST import QueryFeatureMLST
from storage.variant.SampleSet import SampleSet
from storage.api.SamplesQuery import SamplesQuery
from storage.connector.DataIndexConnection import DataIndexConnection


class SamplesQueryIndex(SamplesQuery):

    def __init__(self, connection: DataIndexConnection, sample_set: SampleSet = SamplesQuery.ALL_SAMPLES):
        super().__init__()
        self._query_connection = connection
        self._sample_set = sample_set

    @property
    def sample_set(self) -> SampleSet:
        return self._sample_set

    def _intersect_sample_set(self, other: SampleSet) -> SampleSet:
        if self.sample_set == self.ALL_SAMPLES:
            return other
        elif other == self.ALL_SAMPLES:
            return self.sample_set
        else:
            return self.sample_set.intersection(other)

    def dataframe(self) -> pd.DataFrame:
        pass

    def and_(self, other):
        if isinstance(other, SamplesQuery):
            intersect_set = self._intersect_sample_set(other.sample_set)
            return self._create_from(self._query_connection, intersect_set)
        else:
            raise Exception(f'Cannot perform an "and" on object {other}')

    def __and__(self, other):
        return self.and_(other)

    def __len__(self):
        if self.sample_set == self.ALL_SAMPLES:
            raise Exception('Getting length of all samples is currently unsupported')
        else:
            return len(self.sample_set)

    def has(self, feature: QueryFeature) -> SamplesQuery:
        found_set = self._query_connection.sample_service.find_sample_sets_by_features([feature])
        intersect_found = self._intersect_sample_set(found_set)

        return self._create_from(self._query_connection, intersect_found)

    def has_mutation(self, mutation_feature: str) -> SamplesQuery:
        return self.has(QueryFeatureMutation(mutation_feature))

    def has_allele(self, allele_feature: str) -> SamplesQuery:
        return self.has(QueryFeatureMLST(allele_feature))

    def within(self, distance: float, sample_name: str, within_type: str) -> SamplesQuery:
        raise Exception('Not implemented')

    def within_kmer(self, distance: float, sample_name: str) -> SamplesQuery:
        return self.within(distance, sample_name, 'kmer')

    def is_type(self, sample_type) -> SamplesQuery:
        raise Exception('Not implemented')

    def _create_from(self, connection: DataIndexConnection, sample_set: SampleSet) -> SamplesQuery:
        return SamplesQueryIndex(connection=self._query_connection,
                                 sample_set=sample_set)
