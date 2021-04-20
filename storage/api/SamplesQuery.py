from __future__ import annotations

import abc
from typing import Union

import pandas as pd

from storage.variant.SampleSet import SampleSet
from storage.variant.model.QueryFeature import QueryFeature


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
    def toframe(self) -> pd.DataFrame:
        pass

    @abc.abstractmethod
    def intersect(self, sample_set: SampleSet, query_message: str = None) -> SamplesQuery:
        pass

    @abc.abstractmethod
    def and_(self, other: SamplesQuery) -> SamplesQuery:
        pass

    @abc.abstractmethod
    def build_tree(self, kind: str, scope: str, **kwargs) -> SamplesQuery:
        pass

    @abc.abstractmethod
    def has(self, feature: Union[QueryFeature, str], kind = None) -> SamplesQuery:
        pass

    @abc.abstractmethod
    def within(self, distance: float, sample_name: str, units: str) -> SamplesQuery:
        pass

    @abc.abstractmethod
    def is_type(self, sample_type) -> SamplesQuery:
        pass

    @abc.abstractmethod
    def is_empty(self):
        pass

    @property
    @abc.abstractmethod
    def tree(self):
        pass

    @abc.abstractmethod
    def __and__(self, other):
        pass

    @abc.abstractmethod
    def __len__(self):
        pass
