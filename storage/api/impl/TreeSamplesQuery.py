from __future__ import annotations
import copy

from typing import Union

import pandas as pd
from ete3 import Tree, TreeStyle

from storage.api.SamplesQuery import SamplesQuery
from storage.api.impl.TreeBuilderReferenceMutations import TreeBuilderReferenceMutations
from storage.api.viewer.TreeStyler import TreeStyler, DEFAULT_HIGHLIGHT_STYLES
from storage.connector.DataIndexConnection import DataIndexConnection
from storage.variant.SampleSet import SampleSet
from storage.variant.model.QueryFeature import QueryFeature


class TreeSamplesQuery(SamplesQuery):
    BUILD_TREE_KINDS = ['mutation']
    DISTANCE_UNITS = ['substitutions', 'substitutions/site']

    def __init__(self, connection: DataIndexConnection, wrapped_query: SamplesQuery, tree: Tree,
                 alignment_length: int):
        super().__init__()
        self._query_connection = connection
        self._wrapped_query = wrapped_query
        self._tree = tree
        self._alignment_length = alignment_length

    @property
    def universe_set(self) -> SampleSet:
        return self._wrapped_query.universe_set

    @property
    def sample_set(self) -> SampleSet:
        return self._wrapped_query.sample_set

    def intersect(self, sample_set: SampleSet, query_message: str = None) -> SamplesQuery:
        intersected_query = self._wrapped_query.intersect(sample_set=sample_set, query_message=query_message)
        return self._wrap_create(intersected_query)

    def _wrap_create(self, wrapped_query: SamplesQuery) -> TreeSamplesQuery:
        return TreeSamplesQuery(connection=self._query_connection,
                                wrapped_query=wrapped_query,
                                tree=self._tree,
                                alignment_length=self._alignment_length)

    def toframe(self, exclude_absent: bool = True) -> pd.DataFrame:
        return self._wrapped_query.toframe(exclude_absent=exclude_absent)

    def summary(self) -> pd.DataFrame:
        return self._wrapped_query.summary()

    def and_(self, other: SamplesQuery) -> SamplesQuery:
        return self._wrap_create(self._wrapped_query.and_(other))

    def build_tree(self, kind: str, scope: str, **kwargs):
        return self._wrap_create(self._wrapped_query.build_tree(kind=kind, scope=scope, **kwargs))

    def has(self, feature: Union[QueryFeature, str], kind=None) -> SamplesQuery:
        return self._wrap_create(self._wrapped_query.has(feature=feature, kind=kind))

    def complement(self):
        return self._wrap_create(self._wrapped_query.complement())

    def query_expression(self) -> str:
        return self._wrapped_query.query_expression()

    def tolist(self, names=True):
        return self._wrapped_query.tolist(names=names)

    def within(self, distance: float, sample_name: str, units: str) -> SamplesQuery:
        if units == 'substitutions':
            distance_multiplier = self._alignment_length
        elif units == 'substitutions/site':
            distance_multiplier = 1
        else:
            raise Exception(f'Invalid units=[{units}]. Must be one of {self.DISTANCE_UNITS}')

        sample_service = self._query_connection.sample_service
        sample_name_ids = {s.name: s.id for s in sample_service.find_samples_by_ids(self._wrapped_query.sample_set)}

        sample_leaves = self._tree.get_leaves_by_name(sample_name)
        if len(sample_leaves) != 1:
            raise Exception(
                f'Invalid number of matching leaves for sample [{sample_name}], leaves {sample_leaves}')

        sample_node = sample_leaves[0]

        found_samples_set = set()
        leaves = self._tree.get_leaves()
        for leaf in leaves:
            if leaf.name not in sample_name_ids:
                continue
            sample_distance_to_other_sample = sample_node.get_distance(leaf) * distance_multiplier

            if sample_distance_to_other_sample <= distance:
                found_samples_set.add(sample_name_ids[leaf.name])

        found_samples = SampleSet(found_samples_set)
        return self.intersect(found_samples, f'within({distance} {units} of {sample_name})')

    def tree_styler(self, initial_style: TreeStyle = TreeStyle()) -> TreeStyler:
        return TreeStyler(tree=copy.deepcopy(self._tree),
                          default_highlight_styles=DEFAULT_HIGHLIGHT_STYLES,
                          tree_style=initial_style)

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
        universe_length = len(self.universe_set)
        percent_selected = (len(self) / universe_length) * 100
        return f'<{self.__class__.__name__}[{percent_selected:0.0f}% ({len(self)}/{universe_length}) samples]>'

    def __repr__(self) -> str:
        return str(self)

    @classmethod
    def create(cls, kind: str, scope: str, database_connection: DataIndexConnection,
               wrapped_query: SamplesQuery, **kwargs) -> TreeSamplesQuery:
        if kind == 'mutation':
            tree_builder = TreeBuilderReferenceMutations(database_connection,
                                                         reference_name=scope)
            tree, alignment_length, tree_samples_set = tree_builder.build(wrapped_query.sample_set, **kwargs)

            wrapped_query_tree_set = wrapped_query.intersect(sample_set=tree_samples_set,
                                                             query_message=f'mutation_tree({scope})')
            tree_samples_query = TreeSamplesQuery(connection=database_connection,
                                                  wrapped_query=wrapped_query_tree_set,
                                                  tree=tree,
                                                  alignment_length=alignment_length)
            return tree_samples_query
        else:
            raise Exception(f'Got kind=[{kind}], only the following kinds are supported: {cls.BUILD_TREE_KINDS}')
