from __future__ import annotations
import logging
import abc

import pandas as pd

from storage.variant.service.MutationQueryService import MutationQueryService
from storage.variant.service.QueryService import QueryService
from storage.variant.model.QueryFeature import QueryFeature
from storage.variant.model.QueryFeatureMutation import QueryFeatureMutation
from storage.variant.model.QueryFeatureMLST import QueryFeatureMLST
from storage.variant.SampleSet import SampleSet
from storage.api.SamplesQuery import SamplesQuery


logger = logging.getLogger(__name__)


def query(connection: QueryConnection) -> SamplesQuery:
    return SamplesQuery(connection=connection)








class SamplesQueryIndex(SamplesQuery):
    ALL_SAMPLES = None

    def __init__(self, connection: QueryConnection, sample_set: SampleSet = ALL_SAMPLES):
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

    def _and(self, other):
        if isinstance(other, SamplesQuery):
            intersect_set = self._intersect_sample_set(other.sample_set)
            return SamplesQueryIndex(connection=self._query_connection,
                                       sample_set=intersect_set)
        else:
            raise Exception(f'Cannot perform an "and" on object {other}')

    def __and__(self, other):
        return self._and(other)

    def __len__(self):
        if self.sample_set == self.ALL_SAMPLES:
            raise Exception('Getting length of all samples is currently unsupported')
        else:
            return len(self.sample_set)

    def has(self, feature: QueryFeature) -> SamplesQuery:
        return self._has(feature)

    def has_mutation(self, mutation_feature: str) -> SamplesQuery:
        return self._has(QueryFeatureMutation(mutation_feature))

    def has_allele(self, allele_feature: str) -> SamplesQuery:
        return self._has(QueryFeatureMLST(allele_feature))

    def within(self, distance: float, sample_name: str, within_type: str) -> SamplesQuery:
        raise Exception('Not implemented')

    def within_kmer(self, distance: float, sample_name: str) -> SamplesQuery:
        return self.within(distance, sample_name, 'kmer')

    def is_type(self, sample_type) -> SamplesQuery:
        raise Exception('Not implemented')

    def _has(self, feature: QueryFeature) -> SamplesQuery:
        query_service = self._query_connection.get_feature_service(feature)
        found_set = query_service.find_sample_set_by_feature(feature)
        intersect_found = self._intersect_sample_set(found_set)

        return SamplesQueryIndex(connection=self._query_connection,
                                   sample_set=intersect_found)


# class DataFrameSamplesQuery(SamplesQuery):
#
#     def __init__(self):
