from __future__ import annotations

import abc

import pandas as pd

from storage.variant.SampleSet import SampleSet
from storage.variant.model.QueryFeature import QueryFeature


class SamplesQuery(abc.ABC):
    ALL_SAMPLES = SampleSet.create_all()

    def __init__(self):
        pass

    @property
    @abc.abstractmethod
    def sample_set(self) -> SampleSet:
        pass

    @abc.abstractmethod
    def toframe(self) -> pd.DataFrame:
        pass

    @abc.abstractmethod
    def and_(self, other: SamplesQuery) -> SamplesQuery:
        pass

    @abc.abstractmethod
    def has(self, feature: QueryFeature) -> SamplesQuery:
        pass

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
    def is_empty(self):
        pass

    @abc.abstractmethod
    def __and__(self, other):
        pass

    @abc.abstractmethod
    def __len__(self):
        pass
