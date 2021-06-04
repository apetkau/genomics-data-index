from __future__ import annotations

import logging
from typing import cast, List, Optional

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

    def _join_tree_mutations(self, tree: Tree, kind: str, database_connection: DataIndexConnection,
                             wrapped_query: SamplesQuery,
                             leaf_names: List[str],
                             alignment_length: int,
                             reference_name: str = None) -> TreeSamplesQuery:
        if alignment_length is None:
            raise Exception(f'Must set alignment_length for mutation tree.')
        else:
            include_reference = reference_name is not None
            samples_set = wrapped_query.sample_set
            expected_leaves_number = len(samples_set)

            if include_reference:
                expected_leaves_number += 1

            if len(leaf_names) != expected_leaves_number:
                logger.warning(
                    f'Passed tree has {len(leaf_names)} leaves, but only {expected_leaves_number} match samples in the '
                    f'system. Pruning tree to match samples in system.')
                sample_name_ids = database_connection.sample_service.find_sample_name_ids(leaf_names)
                sample_names = list(sample_name_ids.keys())
                if include_reference:
                    sample_names = sample_names + [reference_name]
                tree.prune(sample_names)

            return self._create_tree_samples_query_from_tree(kind=kind,
                                                             connection=database_connection,
                                                             wrapped_query=wrapped_query,
                                                             tree=tree,
                                                             alignment_length=alignment_length,
                                                             reference_name=reference_name,
                                                             include_reference=include_reference)

    def _join_tree_kmers(self, tree: Tree, kind: str, database_connection: DataIndexConnection,
                         wrapped_query: SamplesQuery,
                         leaf_names: List[str]) -> TreeSamplesQuery:
        samples_set = wrapped_query.sample_set
        expected_leaves_number = len(samples_set)
        if len(leaf_names) != expected_leaves_number:
            logger.warning(
                f'Passed tree has {len(leaf_names)} leaves, but only {expected_leaves_number} match samples in the '
                f'system. Pruning tree to match samples in system.')
            sample_name_ids = database_connection.sample_service.find_sample_name_ids(leaf_names)
            sample_names = list(sample_name_ids.keys())
            tree.prune(sample_names)

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

        leaf_names = tree.get_leaf_names()
        samples_set = SampleSet(database_connection.sample_service.find_sample_name_ids(leaf_names).values())

        wrapped_query_tree_set = wrapped_query.intersect(sample_set=samples_set,
                                                         query_message=f'join_tree({len(leaf_names)} leaves)')

        if kind == 'mutation' or kind == 'mutation_experimental':
            return self._join_tree_mutations(tree=tree, kind=kind, database_connection=database_connection,
                                             wrapped_query=wrapped_query_tree_set,
                                             leaf_names=leaf_names,
                                             **kwargs)
        elif kind == 'kmer' or kind == 'kmers':
            return self._join_tree_kmers(tree=tree, kind=kind, database_connection=database_connection,
                                         wrapped_query=wrapped_query_tree_set,
                                         leaf_names=leaf_names)
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
