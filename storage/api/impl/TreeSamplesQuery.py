from __future__ import annotations
from typing import Union

import pandas as pd
from ete3 import Tree

from storage.api.SamplesQuery import SamplesQuery
from storage.variant.SampleSet import SampleSet
from storage.variant.model.QueryFeature import QueryFeature
from storage.connector.DataIndexConnection import DataIndexConnection
from storage.api.impl.TreeBuilderReferenceMutations import TreeBuilderReferenceMutations


class TreeSamplesQuery(SamplesQuery):

    BUILD_TREE_KINDS = ['mutation']

    def __init__(self, wrapped_query: SamplesQuery, tree: Tree):
        super().__init__()
        self._wrapped_query = wrapped_query
        self._tree = tree

    @property
    def sample_set(self) -> SampleSet:
        return self._wrapped_query.sample_set

    def intersect(self, sample_set: SampleSet, query_message: str = None) -> SamplesQuery:
        intersected_query = self._wrapped_query.intersect(sample_set=sample_set, query_message=query_message)
        return self._wrap_create(intersected_query)

    def _wrap_create(self, wrapped_query: SamplesQuery) -> TreeSamplesQuery:
        return TreeSamplesQuery(wrapped_query=wrapped_query, tree=self._tree)

    def toframe(self) -> pd.DataFrame:
        return self._wrapped_query.toframe()

    def and_(self, other: SamplesQuery) -> SamplesQuery:
        return self._wrap_create(self._wrapped_query.and_(other))

    def build_tree(self, kind: str, scope: str, **kwargs):
        return self._wrap_create(self._wrapped_query.build_tree(kind=kind, scope=scope, **kwargs))

    def has(self, feature: Union[QueryFeature, str], kind=None) -> SamplesQuery:
        return self._wrap_create(self._wrapped_query.has(feature=feature, kind=kind))

    def within(self, distance: float, sample_name: str, kind: str) -> SamplesQuery:
        raise Exception('implement me')

    def is_type(self, sample_type) -> SamplesQuery:
        return self._wrap_create(self._wrapped_query.is_type(sample_type))

    def is_empty(self):
        return self._wrapped_query.is_empty()

    @property
    def tree(self):
        return self._tree

    def __and__(self, other):
        return self.and_(other)

    def __len__(self):
        return len(self._wrapped_query)

    def __str__(self) -> str:
        return f'<TreeSamplesQuery(samples={len(self)})>'

    def __repr__(self) -> str:
        return str(self)

    @classmethod
    def create(cls, kind: str, scope: str, database_connection: DataIndexConnection,
               wrapped_query: SamplesQuery, **kwargs) -> TreeSamplesQuery:
        if kind == 'mutation':
            tree_builder = TreeBuilderReferenceMutations(database_connection,
                                                         reference_name=scope)
            tree, tree_samples_set = tree_builder.build(wrapped_query.sample_set, **kwargs)

            wrapped_query_tree_set = wrapped_query.intersect(sample_set=tree_samples_set,
                                                             query_message=f'mutation_tree({scope})')
            tree_samples_query = TreeSamplesQuery(wrapped_query=wrapped_query_tree_set,
                                                  tree=tree)
            return tree_samples_query
        else:
            raise Exception(f'Got kind=[{kind}], only the following kinds are supported: {cls.BUILD_TREE_KINDS}')
