from __future__ import annotations

import copy
from typing import Union, List, Dict, Set

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
    WITHIN_TYPES = ['distance', 'mrca']

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

    def summary_features(self, kind: str = 'mutations', **kwargs) -> pd.DataFrame:
        return self._wrapped_query.summary_features(kind=kind, **kwargs)

    def tofeaturesset(self, kind: str = 'mutations', selection: str = 'all') -> Set[str]:
        return self._wrapped_query.tofeaturesset(kind=kind, selection=selection)

    def and_(self, other: SamplesQuery) -> SamplesQuery:
        return self._wrap_create(self._wrapped_query.and_(other))

    def build_tree(self, kind: str, scope: str, **kwargs):
        return TreeSamplesQuery.create(kind=kind, scope=scope, database_connection=self._query_connection,
                                       wrapped_query=self, **kwargs)

    def has(self, feature: Union[QueryFeature, str], kind=None) -> SamplesQuery:
        return self._wrap_create(self._wrapped_query.has(feature=feature, kind=kind))

    def complement(self):
        return self._wrap_create(self._wrapped_query.complement())

    def query_expression(self) -> str:
        return self._wrapped_query.query_expression()

    def tolist(self, names=True):
        return self._wrapped_query.tolist(names=names)

    def _within_distance(self, sample_names: Union[str, List[str]], kind: str, distance: float,
                         units: str) -> SamplesQuery:
        if units == 'substitutions':
            distance_multiplier = self._alignment_length
        elif units == 'substitutions/site':
            distance_multiplier = 1
        else:
            raise Exception(f'Invalid units=[{units}]. Must be one of {self.DISTANCE_UNITS}')

        if isinstance(sample_names, list):
            raise NotImplementedError
        elif not isinstance(sample_names, str):
            raise Exception(f'Invalid type for sample_names=[{sample_names}]')

        sample_name_ids = self._get_sample_name_ids()

        sample_leaves = self._tree.get_leaves_by_name(sample_names)
        if len(sample_leaves) != 1:
            raise Exception(
                f'Invalid number of matching leaves for sample [{sample_names}], leaves {sample_leaves}')

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
        return self.intersect(found_samples, f'within({distance} {units} of {sample_names})')

    def _get_sample_name_ids(self) -> Dict[str, int]:
        sample_service = self._query_connection.sample_service
        return {s.name: s.id for s in sample_service.find_samples_by_ids(self._wrapped_query.sample_set)}

    def _within_mrca(self, sample_names: Union[str, List[str]]) -> SamplesQuery:
        if isinstance(sample_names, str):
            sample_names = [sample_names]

        sample_name_ids = self._get_sample_name_ids()
        sample_leaves_list = []
        for name in sample_names:
            sample_leaves = self._tree.get_leaves_by_name(name)
            if len(sample_leaves) != 1:
                raise Exception(
                    f'Invalid number of matching leaves for sample [{name}], leaves {sample_leaves}')
            else:
                sample_leaves_list.append(sample_leaves[0])

        if len(sample_leaves_list) == 0:
            raise Exception(f'Should at least have some leaves in the tree matching sample_names={sample_names}')
        elif len(sample_leaves_list) == 1:
            found_sample_names = sample_names
        else:
            first_sample_leaf = sample_leaves_list.pop()
            ancestor_node = first_sample_leaf.get_common_ancestor(sample_leaves_list)
            found_sample_names = ancestor_node.get_leaf_names()

        found_samples_list = []
        for name in found_sample_names:
            if name in sample_name_ids:
                found_samples_list.append(sample_name_ids[name])
        found_samples = SampleSet(found_samples_list)
        return self.intersect(found_samples, f'within(mrca of {sample_names})')

    def within(self, sample_names: Union[str, List[str]], kind: str = 'distance', **kwargs) -> SamplesQuery:
        if kind == 'distance':
            return self._within_distance(sample_names=sample_names, kind=kind, **kwargs)
        elif kind == 'mrca':
            return self._within_mrca(sample_names=sample_names)
        else:
            raise Exception(f'kind=[{kind}] is not supported. Must be one of {self.WITHIN_TYPES}')

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