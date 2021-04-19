from typing import Union

import pandas as pd

from storage.api.SamplesQuery import SamplesQuery
from storage.variant.SampleSet import SampleSet
from storage.variant.model.QueryFeature import QueryFeature


class TreeSamplesQuery(SamplesQuery):

    def __init__(self, wrapped_query: SamplesQuery):
        super().__init__()
        self._wrapped_query = wrapped_query

    @property
    def sample_set(self) -> SampleSet:
        return self._wrapped_query.sample_set

    def toframe(self) -> pd.DataFrame:
        return self._wrapped_query.toframe()

    def and_(self, other: SamplesQuery) -> SamplesQuery:
        return self._wrapped_query.and_(other)

    def tree(self, kind: str, scope: str, **kwargs):
        return self._wrapped_query.tree(kind=kind, scope=scope, **kwargs)

    def has(self, feature: Union[QueryFeature, str], kind=None) -> SamplesQuery:
        return self._wrapped_query.has(feature=feature, kind=kind)

    def within(self, distance: float, sample_name: str, kind: str) -> SamplesQuery:
        raise Exception('implement me')

    def is_type(self, sample_type) -> SamplesQuery:
        return self._wrapped_query.is_type(sample_type)

    def is_empty(self):
        return self._wrapped_query.is_empty()

    def __and__(self, other):
        return self._wrapped_query and other

    def __len__(self):
        return len(self._wrapped_query)
