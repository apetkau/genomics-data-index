from __future__ import annotations
import copy

from typing import Callable, Set

from ete3 import Tree, TreeNode

from genomics_data_index.api.query.SamplesQuery import SamplesQuery
from genomics_data_index.api.query.impl.WrappedSamplesQuery import WrappedSamplesQuery
from genomics_data_index.configuration.connector.DataIndexConnection import DataIndexConnection

from genomics_data_index.api.query.impl.TreeSamplesQuery import TreeSamplesQuery
from genomics_data_index.storage.SampleSet import SampleSet


class ExperimentalTreeSamplesQuery(TreeSamplesQuery):
    """
    This class contains experimental selection/query statements with respect to a tree
    that I have started implementing but haven't had time to properly test out yet. Eventually
    these may make it into the regular TreeSamplesQuery class.
    """

    def __init__(self, connection: DataIndexConnection, wrapped_query: SamplesQuery, tree: Tree,
                 alignment_length: int, reference_name: str, reference_included: bool):
        super().__init__(connection=connection, wrapped_query=wrapped_query,
                         tree=tree, alignment_length=alignment_length,
                         reference_name=reference_name,
                         reference_included=reference_included)

    def _tree_copy(self):
        return self._tree.copy(method='deepcopy')

    def set_outgroup(self, sample_name: str) -> SamplesQuery:
        tree = self._tree_copy()
        tree.set_outgroup(sample_name)
        return self._create_from_tree_internal(tree)

    def relabel_samples(self, rename_func: Callable[[str], str]) -> SamplesQuery:
        tree = self._tree_copy()
        for node in tree.traverse("postorder"):
            if node.is_leaf():
                node.name = rename_func(node.name)
        return self._create_from_tree_internal(tree)

    def prune(self, preserve_branch_length: bool = True) -> SamplesQuery:
        """
        Prunes tree down to whatever the current query is.
        """
        tree = self._tree_copy()
        tree.prune(self.tolist(), preserve_branch_length=preserve_branch_length)
        return self._create_from_tree_internal(tree)

    def label_internal_nodes(self) -> SamplesQuery:
        tree = self._tree_copy()
        self._label_internal_nodes_recursive(tree)
        return self._create_from_tree_internal(tree)

    def _label_internal_nodes_recursive(self, node: TreeNode) -> Set[str]:
        if node.is_leaf():
            return {node.name}
        else:
            children_names = set()
            children = node.get_children()
            for child in children:
                children_names.update(self._label_internal_nodes_recursive(child))

            node.name = f'{children_names}'
            return children_names

    def _wrap_create(self, wrapped_query: SamplesQuery, universe_set: SampleSet = None) -> WrappedSamplesQuery:
        tree = self._tree_copy()
        tree.prune(wrapped_query.tolist(), preserve_branch_length=True)
        return ExperimentalTreeSamplesQuery(connection=self._query_connection,
                                            wrapped_query=wrapped_query,
                                            tree=tree,
                                            alignment_length=self._alignment_length,
                                            reference_name=self.reference_name,
                                            reference_included=self.reference_included)

    def _create_from_tree_internal(self, tree: Tree) -> SamplesQuery:
        return ExperimentalTreeSamplesQuery(connection=self._query_connection,
                                            wrapped_query=self._wrapped_query,
                                            tree=tree,
                                            alignment_length=self._alignment_length,
                                            reference_name=self.reference_name,
                                            reference_included=self.reference_included)

    @classmethod
    def from_tree_query(self, query: TreeSamplesQuery) -> ExperimentalTreeSamplesQuery:
        """
        Creates new ExperimentalTreeSamplesQuery to match the same information as the passed query.
        :param query: The query to copy from.
        :return: An equivalent experimental query.
        """
        return ExperimentalTreeSamplesQuery(connection=query._query_connection,
                                            wrapped_query=query._wrapped_query,
                                            tree=query.tree,
                                            alignment_length=query._alignment_length,
                                            reference_name=query.reference_name,
                                            reference_included=query.reference_included)
