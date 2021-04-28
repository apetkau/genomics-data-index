from __future__ import annotations
import copy

from typing import Callable

from ete3 import Tree

from genomics_data_index.api.query.SamplesQuery import SamplesQuery
from genomics_data_index.configuration.connector.DataIndexConnection import DataIndexConnection

from genomics_data_index.api.query.impl.TreeSamplesQuery import TreeSamplesQuery


class ExperimentalTreeSamplesQuery(TreeSamplesQuery):
    """
    This class contains experimental selection/query statements with respect to a tree
    that I have started implementing but haven't had time to properly test out yet. Eventually
    these may make it into the regular TreeSamplesQuery class.
    """

    def __init__(self, connection: DataIndexConnection, wrapped_query: SamplesQuery, tree: Tree,
                 alignment_length: int):
        super().__init__(connection=connection, wrapped_query=wrapped_query,
                         tree=tree, alignment_length=alignment_length)

    def relabel_samples(self, rename_func: Callable[[str], str]) -> SamplesQuery:
        tree = copy.deepcopy(self._tree)
        for node in tree.traverse("postorder"):
            if node.is_leaf():
                node.name = rename_func(node.name)
        return self._create_from_tree_internal(tree)

    def _create_from_tree_internal(self, tree: Tree) -> SamplesQuery:
        return ExperimentalTreeSamplesQuery(connection=self._query_connection,
                                            wrapped_query=self._wrapped_query,
                                            tree=tree,
                                            alignment_length=self._alignment_length)


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
                                            alignment_length=query._alignment_length)
