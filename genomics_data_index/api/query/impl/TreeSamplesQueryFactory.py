from __future__ import annotations

from genomics_data_index.api.query.SamplesQuery import SamplesQuery
from genomics_data_index.api.query.impl.TreeBuilderReferenceMutations import TreeBuilderReferenceMutations
from genomics_data_index.api.query.impl.TreeSamplesQuery import TreeSamplesQuery
from genomics_data_index.configuration.connector.DataIndexConnection import DataIndexConnection


class TreeSamplesQueryFactory:
    """
    Builds TreeSamplesQuery objects.
    Used mainly to switch between regular and experimental tree queries while avoiding circular dependencies.
    """
    BUILD_TREE_KINDS = ['mutation']

    _instance = None

    def __init__(self):
        pass

    def build_tree(self, kind: str, scope: str, database_connection: DataIndexConnection,
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
            raise Exception(f'Got kind=[{kind}], only the following kinds are supported: {self.BUILD_TREE_KINDS}')

    @classmethod
    def instance(cls) -> TreeSamplesQueryFactory:
        if cls._instance is None:
            cls._instance = TreeSamplesQueryFactory()
        return cls._instance
