from typing import Union, List

from ete3 import Tree

from genomics_data_index.api.query.SamplesQuery import SamplesQuery
from genomics_data_index.api.query.impl.TreeBuilderReferenceMutations import TreeBuilderReferenceMutations
from genomics_data_index.api.query.impl.WrappedSamplesQuery import WrappedSamplesQuery
from genomics_data_index.configuration.connector.DataIndexConnection import DataIndexConnection

from genomics_data_index.api.query.impl.TreeSamplesQuery import TreeSamplesQuery
from genomics_data_index.storage.SampleSet import SampleSet


class MutationTreeSamplesQuery(TreeSamplesQuery):

    DISTANCE_UNITS = ['substitutions', 'substitutions/site']

    def __init__(self, connection: DataIndexConnection, wrapped_query: SamplesQuery, tree: Tree,
                 alignment_length: int, reference_name: str, reference_included: bool):
        super().__init__(connection=connection, wrapped_query=wrapped_query, tree=tree)
        self._alignment_length = alignment_length
        self._reference_name = reference_name
        self._reference_included = reference_included

    @property
    def reference_name(self):
        return self._reference_name

    @property
    def reference_included(self):
        return self._reference_included

    def _wrap_create(self, wrapped_query: SamplesQuery, universe_set: SampleSet = None) -> WrappedSamplesQuery:
        return MutationTreeSamplesQuery(connection=self._query_connection,
                                        wrapped_query=wrapped_query,
                                        tree=self._tree,
                                        alignment_length=self._alignment_length,
                                        reference_name=self._reference_name,
                                        reference_included=self._reference_included)

    def _within_distance_internal(self, sample_names: Union[str, List[str]], distance: float,
                         units: str) -> SamplesQuery:
        if units == 'substitutions':
            distance_multiplier = self._alignment_length
        elif units == 'substitutions/site':
            distance_multiplier = 1
        else:
            raise Exception(f'Invalid units=[{units}]. Must be one of {self.DISTANCE_UNITS}')

        if isinstance(sample_names, list):
            raise NotImplementedError
        elif not isinstance(sample_names, str):
            raise Exception(f'Invalid type for sample_names=[{sample_names}]')

        sample_name_ids = self._get_sample_name_ids()

        sample_leaves = self._tree.get_leaves_by_name(sample_names)
        if len(sample_leaves) != 1:
            raise Exception(
                f'Invalid number of matching leaves for sample [{sample_names}], leaves {sample_leaves}')

        sample_node = sample_leaves[0]

        found_samples_set = set()
        leaves = self._tree.get_leaves()
        for leaf in leaves:
            if leaf.name not in sample_name_ids:
                continue
            sample_distance_to_other_sample = sample_node.get_distance(leaf) * distance_multiplier

            if sample_distance_to_other_sample <= distance:
                found_samples_set.add(sample_name_ids[leaf.name])

        found_samples = SampleSet(found_samples_set)
        return self.intersect(found_samples, f'within({distance} {units} of {sample_names})')

    def _can_handle_distance_units(self, units: str) -> bool:
        return units in self.DISTANCE_UNITS

    def _distance_units(self) -> List[str]:
        return self.DISTANCE_UNITS

    def build_tree(self, kind: str, scope: str, **kwargs):
        return MutationTreeSamplesQuery.create(kind=kind, scope=scope, database_connection=self._query_connection,
                                               wrapped_query=self, **kwargs)

    @classmethod
    def create(cls, kind: str, scope: str, database_connection: DataIndexConnection,
               wrapped_query: SamplesQuery, include_reference=True, **kwargs) -> TreeSamplesQuery:
        if kind == 'mutation':
            tree_builder = TreeBuilderReferenceMutations(database_connection,
                                                         reference_name=scope)
            tree, alignment_length, tree_samples_set = tree_builder.build(wrapped_query.sample_set,
                                                                          include_reference=include_reference,
                                                                          **kwargs)

            wrapped_query_tree_set = wrapped_query.intersect(sample_set=tree_samples_set,
                                                             query_message=f'mutation_tree({scope})')
            tree_samples_query = MutationTreeSamplesQuery(connection=database_connection,
                                                          wrapped_query=wrapped_query_tree_set,
                                                          tree=tree,
                                                          alignment_length=alignment_length,
                                                          reference_name=scope,
                                                          reference_included=include_reference)
            return tree_samples_query
        else:
            raise Exception(f'Got kind=[{kind}], only the following kinds are supported: {cls.BUILD_TREE_KINDS}')
