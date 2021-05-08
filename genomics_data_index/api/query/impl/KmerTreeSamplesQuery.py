from typing import cast, List, Union

from ete3 import ClusterTree

from genomics_data_index.api.query.SamplesQuery import SamplesQuery
from genomics_data_index.api.query.impl.TreeSamplesQuery import TreeSamplesQuery
from genomics_data_index.api.query.impl.WrappedSamplesQuery import WrappedSamplesQuery
from genomics_data_index.configuration.connector.DataIndexConnection import DataIndexConnection
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

    # Override primarily to specify a default value for 'units'
    def _within_distance(self, sample_names: Union[str, List[str]], distance: float,
                         units: str = 'kmer_jaccard', **kwargs) -> SamplesQuery:
        return super()._within_distance(sample_names=sample_names, distance=distance,
                                        units=units, **kwargs)

    def _within_distance_internal(self, sample_names: Union[str, List[str]], distance: float,
                                  units: str) -> SamplesQuery:
        raise NotImplementedError(f'Not implemented for {KmerTreeSamplesQuery.__class__}')

    def _can_handle_distance_units(self, units: str) -> bool:
        return False

    def _distance_units(self) -> List[str]:
        return self._wrapped_query._distance_units()

    def build_tree(self, kind: str, scope: str, **kwargs):
        raise Exception('TODO: implement this')
