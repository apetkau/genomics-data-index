from __future__ import annotations

import copy
from typing import Union, List

from ete3 import Tree, TreeStyle

from storage.api.SamplesQuery import SamplesQuery
from storage.api.impl.WrappedSamplesQuery import WrappedSamplesQuery
from storage.api.impl.TreeBuilderReferenceMutations import TreeBuilderReferenceMutations
from storage.api.viewer.TreeStyler import TreeStyler, DEFAULT_HIGHLIGHT_STYLES
from storage.configuration.connector import DataIndexConnection
from storage.variant.SampleSet import SampleSet


class TreeSamplesQuery(WrappedSamplesQuery):
    BUILD_TREE_KINDS = ['mutation']
    DISTANCE_UNITS = ['substitutions', 'substitutions/site']
    WITHIN_TYPES = ['distance', 'mrca']

    def __init__(self, connection: DataIndexConnection, wrapped_query: SamplesQuery, tree: Tree,
                 alignment_length: int):
        super().__init__(connection=connection, wrapped_query=wrapped_query)
        self._tree = tree
        self._alignment_length = alignment_length

    def _wrap_create(self, wrapped_query: SamplesQuery) -> WrappedSamplesQuery:
        return TreeSamplesQuery(connection=self._query_connection,
                                wrapped_query=wrapped_query,
                                tree=self._tree,
                                alignment_length=self._alignment_length)

    def build_tree(self, kind: str, scope: str, **kwargs):
        return TreeSamplesQuery.create(kind=kind, scope=scope, database_connection=self._query_connection,
                                       wrapped_query=self, **kwargs)

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

    @property
    def tree(self):
        return self._tree

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
