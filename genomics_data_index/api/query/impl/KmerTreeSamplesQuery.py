from ete3 import ClusterTree

from genomics_data_index.api.query.SamplesQuery import SamplesQuery
from genomics_data_index.configuration.connector.DataIndexConnection import DataIndexConnection

from genomics_data_index.api.query.impl.TreeSamplesQuery import TreeSamplesQuery


class KmerTreeSamplesQuery(TreeSamplesQuery):

    def __init__(self, connection: DataIndexConnection, wrapped_query: SamplesQuery, tree: ClusterTree):
        super().__init__(connection=connection, wrapped_query=wrapped_query,
                         tree=tree, alignment_length=None,
                         reference_name=None,
                         reference_included=None)
