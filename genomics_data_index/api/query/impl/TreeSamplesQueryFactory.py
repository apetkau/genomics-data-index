from __future__ import annotations

import logging
from typing import cast, Optional, Set

from ete3 import Tree, ClusterTree

from genomics_data_index.api.query.SamplesQuery import SamplesQuery
from genomics_data_index.api.query.impl.ExperimentalTreeSamplesQuery import ExperimentalTreeSamplesQuery
from genomics_data_index.api.query.impl.KmerTreeSamplesQuery import KmerTreeSamplesQuery
from genomics_data_index.api.query.impl.MutationTreeSamplesQuery import MutationTreeSamplesQuery
from genomics_data_index.api.query.impl.TreeBuilderKmers import TreeBuilderKmers
from genomics_data_index.api.query.impl.TreeBuilderReferenceMutations import TreeBuilderReferenceMutations
from genomics_data_index.api.query.impl.TreeSamplesQuery import TreeSamplesQuery
from genomics_data_index.configuration.connector.DataIndexConnection import DataIndexConnection
from genomics_data_index.storage.SampleSet import SampleSet
from genomics_data_index.storage.model.db import Reference

logger = logging.getLogger(__name__)


class TreeSamplesQueryFactory:
    """
    Builds TreeSamplesQuery objects.
    Used mainly to switch between regular and experimental tree queries while avoiding circular dependencies.
    """
    BUILD_TREE_KINDS = ['mutation', 'mutations', 'mutation_experimental', 'mutations_experimental',
                        'kmer', 'kmers']

    _instance = None

    def __init__(self):
        pass

    def _prune_tree_to_query_samples(self, tree: Tree,
                                     leaf_names: Set[str],
                                     leaf_names_in_query: Set[str],
                                     reference_name: Optional[str]) -> Tree:
        include_reference = reference_name is not None
        leaves_missing_in_database = leaf_names - leaf_names_in_query

        if len(leaves_missing_in_database) == 0:
            if include_reference:
                logger.warning(f'Reference genome name [{reference_name}] is the same as a sample name already '
                               f'in the database.')
        else:
            if not include_reference or leaves_missing_in_database != {reference_name}:
                if include_reference:
                    leaves_to_prune = leaf_names_in_query.copy()
                    leaves_to_prune.add(reference_name)

                    leaves_missing_not_reference = leaves_missing_in_database.copy()

                    if reference_name in leaves_missing_not_reference:
                        leaves_missing_not_reference.remove(reference_name)

                    extra_msg = ' (and not the reference genome name)'
                else:
                    leaves_to_prune = leaf_names_in_query
                    leaves_missing_not_reference = leaves_missing_in_database

                    extra_msg = ''

                if len(leaves_missing_not_reference) > 5:
                    missing_str = ', '.join(list(leaves_missing_not_reference)[0:5])
                    missing_str += ', ...'
                else:
                    missing_str = ', '.join(list(leaves_missing_not_reference))

                logger.warning(
                    f'{len(leaves_missing_not_reference)}/{len(leaf_names)} leaves in the tree are'
                    f' not found in the query{extra_msg}: [{missing_str}].'
                    f' Pruning tree to contain only those samples in the query ('
                    f'{len(leaves_to_prune)}/{len(leaf_names)}).')

                tree = tree.copy(method='newick')
                tree.prune(leaves_to_prune, preserve_branch_length=True)

        return tree

    def _join_tree_mutations(self, tree: Tree, kind: str, database_connection: DataIndexConnection,
                             wrapped_query: SamplesQuery,
                             leaf_names: Set[str],
                             leaf_names_in_query: Set[str],
                             alignment_length: int,
                             reference_name: str = None) -> TreeSamplesQuery:
        if alignment_length is None:
            raise Exception(f'Must set alignment_length for mutation tree.')
        else:
            tree = self._prune_tree_to_query_samples(tree=tree,
                                                     leaf_names=leaf_names,
                                                     leaf_names_in_query=leaf_names_in_query,
                                                     reference_name=reference_name)

            return self._create_tree_samples_query_from_tree(kind=kind,
                                                             connection=database_connection,
                                                             wrapped_query=wrapped_query,
                                                             tree=tree,
                                                             alignment_length=alignment_length,
                                                             reference_name=reference_name,
                                                             include_reference=reference_name is not None)

    def _join_tree_kmers(self, tree: Tree, kind: str, database_connection: DataIndexConnection,
                         wrapped_query: SamplesQuery,
                         leaf_names: Set[str],
                         leaf_names_in_query: Set[str]) -> TreeSamplesQuery:
        tree = self._prune_tree_to_query_samples(tree=tree,
                                                 leaf_names=leaf_names,
                                                 leaf_names_in_query=leaf_names_in_query,
                                                 reference_name=None)

        if not isinstance(tree, ClusterTree):
            logger.warning(f'Passed tree={tree} is not a {ClusterTree.__name__}. '
                           f'Converting to a {ClusterTree.__name__}')
            newick = tree.write()
            cluster_tree = ClusterTree(newick=newick)
        else:
            cluster_tree = cast(ClusterTree, tree)

        return self._create_tree_samples_query_from_tree(kind=kind,
                                                         connection=database_connection,
                                                         wrapped_query=wrapped_query,
                                                         tree=cluster_tree,
                                                         alignment_length=None,
                                                         reference_name=None,
                                                         include_reference=False)

    def join_tree(self, tree: Tree, kind: str, database_connection: DataIndexConnection,
                  wrapped_query: SamplesQuery, **kwargs) -> TreeSamplesQuery:
        if kind not in self.BUILD_TREE_KINDS:
            raise Exception(f'Invalid kind=[{kind}]. Must be one of {self.BUILD_TREE_KINDS}.')

        tree = tree.copy(method='newick')

        leaf_names = set(tree.get_leaf_names())
        sample_name_ids = database_connection.sample_service.find_sample_name_ids(leaf_names)
        samples_set = SampleSet(sample_name_ids.values())

        wrapped_query_tree_set = wrapped_query.intersect(sample_set=samples_set,
                                                         query_message=f'join_tree({len(leaf_names)} leaves)')

        leaf_names_in_query = wrapped_query_tree_set.toset(names=True)

        if kind == 'mutation' or kind == 'mutation_experimental':
            return self._join_tree_mutations(tree=tree, kind=kind, database_connection=database_connection,
                                             wrapped_query=wrapped_query_tree_set,
                                             leaf_names=leaf_names,
                                             leaf_names_in_query=leaf_names_in_query,
                                             **kwargs)
        elif kind == 'kmer' or kind == 'kmers':
            return self._join_tree_kmers(tree=tree, kind=kind, database_connection=database_connection,
                                         wrapped_query=wrapped_query_tree_set,
                                         leaf_names=leaf_names,
                                         leaf_names_in_query=leaf_names_in_query)
        else:
            raise Exception(f'Invalid kind=[{kind}]. Must be one of {self.BUILD_TREE_KINDS}.')

    def build_tree(self, kind: str, database_connection: DataIndexConnection,
                   wrapped_query: SamplesQuery, scope: str = None, include_reference: bool = True,
                   **kwargs) -> TreeSamplesQuery:
        if kind == 'mutation' or kind == 'mutation_experimental':
            if scope is None:
                raise Exception(f'Invalid scope=[{scope}]. You must specify a scope (i.e., reference genome name) for '
                                f'kind=[{kind}].')

            tree_builder = TreeBuilderReferenceMutations(database_connection,
                                                         reference_name=scope)
            tree, alignment_length, tree_samples_set = tree_builder.build(wrapped_query.sample_set,
                                                                          include_reference=include_reference,
                                                                          **kwargs)
        elif kind == 'kmer' or kind == 'kmers':
            tree_builder = TreeBuilderKmers(database_connection)
            tree, alignment_length, tree_samples_set = tree_builder.build(wrapped_query.sample_set,
                                                                          **kwargs)
        else:
            raise Exception(f'Got kind=[{kind}], only the following kinds are supported: {self.BUILD_TREE_KINDS}')

        wrapped_query_tree_set = wrapped_query.intersect(sample_set=tree_samples_set,
                                                         query_message=f'mutation_tree({scope})')

        return self._create_tree_samples_query_from_tree(kind=kind,
                                                         connection=database_connection,
                                                         wrapped_query=wrapped_query_tree_set,
                                                         tree=tree,
                                                         alignment_length=alignment_length,
                                                         reference_name=scope,
                                                         include_reference=include_reference)

    def _create_tree_samples_query_from_tree(self, kind: str, tree: Tree,
                                             alignment_length: Optional[int],
                                             reference_name: Optional[str],
                                             connection: DataIndexConnection,
                                             wrapped_query: SamplesQuery,
                                             include_reference: bool = True) -> TreeSamplesQuery:
        if kind == 'mutation' or kind == 'mutations':
            return MutationTreeSamplesQuery(connection=connection,
                                            wrapped_query=wrapped_query,
                                            tree=tree,
                                            alignment_length=alignment_length,
                                            reference_name=reference_name,
                                            reference_included=include_reference)
        elif kind == 'mutation_experimental' or kind == 'mutations_experimental':
            return ExperimentalTreeSamplesQuery(connection=connection,
                                                wrapped_query=wrapped_query,
                                                tree=tree,
                                                universe_tree=tree,
                                                alignment_length=alignment_length,
                                                reference_name=reference_name,
                                                reference_included=include_reference)
        elif kind == 'kmer' or kind == 'kmers':
            if isinstance(tree, ClusterTree):
                tree = cast(ClusterTree, tree)
            else:
                raise Exception(
                    f'Invalid type for tree=[{tree}]. Expected [{ClusterTree.__name__}], got [{type(tree)}]')

            return KmerTreeSamplesQuery(connection=connection,
                                        wrapped_query=wrapped_query,
                                        tree=tree)
        else:
            raise Exception(f'Got kind=[{kind}], only the following kinds are supported: {self.BUILD_TREE_KINDS}')

    def create_from_reference_genome(self, kind: str, reference_genome: Reference,
                                     connection: DataIndexConnection,
                                     wrapped_query: SamplesQuery,
                                     include_reference: bool = True):
        return self._create_tree_samples_query_from_tree(kind=kind,
                                                         connection=connection, wrapped_query=wrapped_query,
                                                         tree=reference_genome.tree,
                                                         alignment_length=reference_genome.tree_alignment_length,
                                                         reference_name=reference_genome.name,
                                                         include_reference=include_reference)

    @classmethod
    def instance(cls) -> TreeSamplesQueryFactory:
        if cls._instance is None:
            cls._instance = TreeSamplesQueryFactory()
        return cls._instance
