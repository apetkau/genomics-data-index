from typing import cast
from ete3 import ClusterTree

from genomics_data_index.api.query.SamplesQuery import SamplesQuery
from genomics_data_index.api.query.impl.WrappedSamplesQuery import WrappedSamplesQuery
from genomics_data_index.configuration.connector.DataIndexConnection import DataIndexConnection

from genomics_data_index.api.query.impl.TreeSamplesQuery import TreeSamplesQuery
from genomics_data_index.storage.SampleSet import SampleSet


class KmerTreeSamplesQuery(TreeSamplesQuery):

    def __init__(self, connection: DataIndexConnection, wrapped_query: SamplesQuery, tree: ClusterTree):
        super().__init__(connection=connection, wrapped_query=wrapped_query, tree=tree)

    def _wrap_create(self, wrapped_query: SamplesQuery, universe_set: SampleSet = None) -> WrappedSamplesQuery:
        tree = self.tree
        if isinstance(tree, ClusterTree):
            tree = cast(ClusterTree, self.tree)
        else:
            raise Exception(f'Incorrect type of tree [{tree}]. Expected [{ClusterTree.__class__}], got [{type(tree)}].')

        return KmerTreeSamplesQuery(connection=self._query_connection,
                                    wrapped_query=wrapped_query,
                                    tree=tree)

    def build_tree(self, kind: str, scope: str, **kwargs):
        raise Exception('TODO: implement this')
