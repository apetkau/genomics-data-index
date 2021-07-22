from __future__ import annotations

from typing import Callable, Set

from ete3 import Tree, TreeNode

from genomics_data_index.api.query.SamplesQuery import SamplesQuery
from genomics_data_index.api.query.impl.MutationTreeSamplesQuery import MutationTreeSamplesQuery
from genomics_data_index.api.query.impl.WrappedSamplesQuery import WrappedSamplesQuery
from genomics_data_index.configuration.connector.DataIndexConnection import DataIndexConnection
from genomics_data_index.storage.SampleSet import SampleSet


class ExperimentalTreeSamplesQuery(MutationTreeSamplesQuery):
    """
    This class contains experimental selection/query statements with respect to a tree
    that I have started implementing but haven't had time to properly test out yet. Eventually
    these may make it into the regular TreeSamplesQuery class.
    """

    def __init__(self, connection: DataIndexConnection, wrapped_query: SamplesQuery,
                 tree: Tree, universe_tree: Tree,
                 alignment_length: int, reference_name: str, reference_included: bool):
        """
        Builds a new ExperimentalTreeSamplesQuery from the given information. In most normal operations SamplesQuery objects
        are not created directly but are instead created from an :py:class:`genomics_data_index.api.GenomicsDataIndex`
        object or from operations applied to a SamplesQuery (e.g., build_tree() or join_tree()).

        :param connection: A connection to a database containing samples.
        :param wrapped_query: The SamplesQuery to wrap around and join to the passed tree.
        :param tree: The tree to join to this query.
        :param alignment_length: The length of the alignment of mutations used to generate the tree
                                 (used to convert substitutions/site distances to substitutions distances).
        :param universe_tree: A tree representing the entire universe of samples in this query.
        :param reference_name: The name of the reference genome in the tree.
        :param reference_included: True if the reference genome is one of the leaves of the tree, False otherwise.
        :return: A new ExperimentalTreeSamplesQuery object.
        """
        super().__init__(connection=connection, wrapped_query=wrapped_query,
                         tree=tree, alignment_length=alignment_length,
                         reference_name=reference_name,
                         reference_included=reference_included)
        self._universe_tree = universe_tree

    @property
    def universe_tree(self) -> Tree:
        return self._universe_tree

    def reset_universe(self, include_unknown: bool = True) -> SamplesQuery:
        wrapped_reset_universe = self._wrapped_query.reset_universe(include_unknown=include_unknown)
        universe_tree = self._tree_copy_prune(from_query=wrapped_reset_universe, preserve_branch_length=True)
        return ExperimentalTreeSamplesQuery(connection=self._query_connection,
                                            wrapped_query=wrapped_reset_universe,
                                            tree=universe_tree,
                                            universe_tree=universe_tree,
                                            alignment_length=self._alignment_length,
                                            reference_name=self.reference_name,
                                            reference_included=self.reference_included)

    def _tree_copy(self) -> Tree:
        return self._tree.copy(method='cpickle')

    def _universe_tree_copy(self) -> Tree:
        return self._universe_tree.copy(method='cpickle')

    def _tree_copy_prune(self, from_query: SamplesQuery = None, preserve_branch_length: bool = True) -> Tree:
        tree = self._universe_tree_copy()
        if from_query is not None:
            query = from_query
        else:
            query = self

        if self.reference_included:
            nodes_to_keep = query.tolist() + [self.reference_name]
        else:
            nodes_to_keep = query.tolist()
        tree.prune(nodes_to_keep, preserve_branch_length=preserve_branch_length)
        return tree

    def set_outgroup(self, sample_name: str) -> SamplesQuery:
        """
        Sets the outgroup of the tree to the specified sample name.
        :param sample_name: The name to set as an outgroup. Must exist as a leaf in the tree.
        :return: A new ExperimentalTreeSamplesQuery with the defined outgroup set.
        """
        tree = self._tree_copy()
        tree.set_outgroup(sample_name)
        return self._create_from_tree_internal(tree)

    def relabel_samples(self, rename_func: Callable[[str], str]) -> SamplesQuery:
        """
        Relabels samples/leaves in the tree based on the passed function.
        :param rename_func: A function used for relabeling. This function takes as input the leaf name
                            and returns a new name for this leaf.
        :return: A ExperimentalTreeSamplesQuery where the tree has been relabeled.
        """
        tree = self._tree_copy()
        for node in tree.traverse("postorder"):
            if node.is_leaf():
                node.name = rename_func(node.name)
        return self._create_from_tree_internal(tree)

    def prune(self, preserve_branch_length: bool = True) -> SamplesQuery:
        """
        Prunes the tree joined to this ExperimentalTreeSamplesQuery down to only the currently selected samples.
        :param preserve_branch_length: True if branch lengths between the subset of samples should be preserved, False otherwise.
        :return: A new ExperimentalTreeSamplesQuery with the tree pruned down to only the selected samples (and the reference
                 genome if it is one of the leaves of the tree).
        """
        tree = self._tree_copy_prune(preserve_branch_length=preserve_branch_length)
        return self._create_from_tree_internal(tree)

    def label_internal_nodes(self) -> SamplesQuery:
        """
        Relabels internal nodes based on the selected samples.
        :return: A new ExperimentalTreeSamplesQuery with the internal nodes of the tree labeled.
        """
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
        tree = self._tree_copy_prune(from_query=wrapped_query, preserve_branch_length=True)
        return ExperimentalTreeSamplesQuery(connection=self._query_connection,
                                            wrapped_query=wrapped_query,
                                            tree=tree,
                                            universe_tree=self._universe_tree,
                                            alignment_length=self._alignment_length,
                                            reference_name=self.reference_name,
                                            reference_included=self.reference_included)

    def _create_from_tree_internal(self, tree: Tree) -> SamplesQuery:
        return ExperimentalTreeSamplesQuery(connection=self._query_connection,
                                            wrapped_query=self._wrapped_query,
                                            tree=tree,
                                            universe_tree=self._universe_tree,
                                            alignment_length=self._alignment_length,
                                            reference_name=self.reference_name,
                                            reference_included=self.reference_included)

    @classmethod
    def from_tree_query(self, query: MutationTreeSamplesQuery) -> ExperimentalTreeSamplesQuery:
        """
        Creates new ExperimentalTreeSamplesQuery to match the same information as the passed query.
        :param query: The query to copy from.
        :return: An equivalent experimental query.
        """
        return ExperimentalTreeSamplesQuery(connection=query._query_connection,
                                            wrapped_query=query._wrapped_query,
                                            tree=query.tree,
                                            universe_tree=query.tree,
                                            alignment_length=query._alignment_length,
                                            reference_name=query.reference_name,
                                            reference_included=query.reference_included)
