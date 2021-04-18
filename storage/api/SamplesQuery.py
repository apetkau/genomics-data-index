from __future__ import annotations
import abc

import pandas as pd

from storage.variant.model.QueryFeature import QueryFeature
from storage.variant.SampleSet import SampleSet


class SamplesQuery(abc.ABC):
    ALL_SAMPLES = None

    def __init__(self):
        pass

    @property
    @abc.abstractmethod
    def sample_set(self) -> SampleSet:
        pass

    @abc.abstractmethod
    def dataframe(self) -> pd.DataFrame:
        pass

    @abc.abstractmethod
    def _and(self, other: SamplesQuery) -> SamplesQuery:
        pass

    @abc.abstractmethod
    def has(self, feature: QueryFeature) -> SamplesQuery:
        return self._has(feature)

    @abc.abstractmethod
    def has_mutation(self, mutation_feature: str) -> SamplesQuery:
        pass

    @abc.abstractmethod
    def has_allele(self, allele_feature: str) -> SamplesQuery:
        pass

    @abc.abstractmethod
    def within(self, distance: float, sample_name: str, within_type: str) -> SamplesQuery:
        pass

    @abc.abstractmethod
    def within_kmer(self, distance: float, sample_name: str) -> SamplesQuery:
        pass

    @abc.abstractmethod
    def is_type(self, sample_type) -> SamplesQuery:
        pass

    @abc.abstractmethod
    def _has(self, feature: QueryFeature) -> SamplesQuery:
        pass

    @abc.abstractmethod
    def __and__(self, other):
        pass

    @abc.abstractmethod
    def __len__(self):
        pass
