from __future__ import annotations

from ete3 import Tree

from genomics_data_index.api.query.SamplesQuery import SamplesQuery
from genomics_data_index.api.query.impl.ExperimentalTreeSamplesQuery import ExperimentalTreeSamplesQuery
from genomics_data_index.api.query.impl.TreeBuilderReferenceMutations import TreeBuilderReferenceMutations
from genomics_data_index.api.query.impl.TreeSamplesQuery import TreeSamplesQuery
from genomics_data_index.configuration.connector.DataIndexConnection import DataIndexConnection
from genomics_data_index.storage.model.db import Reference


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

    def build_tree(self, kind: str, scope: str, database_connection: DataIndexConnection,
                   wrapped_query: SamplesQuery, include_reference: bool = True, **kwargs) -> TreeSamplesQuery:
        if kind == 'mutation' or kind == 'mutation_experimental':
            tree_builder = TreeBuilderReferenceMutations(database_connection,
                                                         reference_name=scope)
            tree, alignment_length, tree_samples_set = tree_builder.build(wrapped_query.sample_set,
                                                                          include_reference=include_reference,
                                                                          **kwargs)

            wrapped_query_tree_set = wrapped_query.intersect(sample_set=tree_samples_set,
                                                             query_message=f'mutation_tree({scope})')
        else:
            raise Exception(f'Got kind=[{kind}], only the following kinds are supported: {self.BUILD_TREE_KINDS}')

        return self._create_tree_samples_query_from_tree(kind=kind,
                                                         connection=database_connection,
                                                         wrapped_query=wrapped_query_tree_set,
                                                         tree=tree,
                                                         alignment_length=alignment_length,
                                                         reference_name=scope,
                                                         include_reference=include_reference)

    def _create_tree_samples_query_from_tree(self, kind: str, tree: Tree,
                                             alignment_length: int,
                                             reference_name: str,
                                             connection: DataIndexConnection,
                                             wrapped_query: SamplesQuery,
                                             include_reference: bool = True) -> TreeSamplesQuery:
        if kind == 'mutation' or kind == 'mutations':
            return TreeSamplesQuery(connection=connection,
                                    wrapped_query=wrapped_query,
                                    tree=tree,
                                    alignment_length=alignment_length,
                                    reference_name=reference_name,
                                    reference_included=include_reference)
        elif kind == 'mutation_experimental' or kind == 'mutations_experimental':
            return ExperimentalTreeSamplesQuery(connection=connection,
                                                wrapped_query=wrapped_query,
                                                tree=tree,
                                                alignment_length=alignment_length,
                                                reference_name=reference_name,
                                                reference_included=include_reference)
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
