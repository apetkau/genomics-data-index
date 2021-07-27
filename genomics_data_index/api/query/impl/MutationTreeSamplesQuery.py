import logging
from typing import Union, List

from ete3 import Tree

from genomics_data_index.api.query.SamplesQuery import SamplesQuery
from genomics_data_index.api.query.impl.TreeSamplesQuery import TreeSamplesQuery
from genomics_data_index.api.query.impl.WrappedSamplesQuery import WrappedSamplesQuery
from genomics_data_index.configuration.connector.DataIndexConnection import DataIndexConnection
from genomics_data_index.storage.SampleSet import SampleSet

logger = logging.getLogger(__name__)


class MutationTreeSamplesQuery(TreeSamplesQuery):
    """
    A MutationTreeSamplesQuery represents a query with a tree derived from nucleotide mutation data.
    """

    DISTANCE_UNITS = ['substitutions', 'substitutions/site']

    def __init__(self, connection: DataIndexConnection, wrapped_query: SamplesQuery, tree: Tree,
                 alignment_length: int, reference_name: str, reference_included: bool):
        """
        Builds a new MutationTreeSamplesQuery from the given information. In most normal operations SamplesQuery objects
        are not created directly but are instead created from an :py:class:`genomics_data_index.api.GenomicsDataIndex`
        object or from operations applied to a SamplesQuery (e.g., build_tree() or join_tree()).

        :param connection: A connection to a database containing samples.
        :param wrapped_query: The SamplesQuery to wrap around and join to the passed tree.
        :param tree: The tree to join to this query.
        :param alignment_length: The length of the alignment of mutations used to generate the tree
                                 (used to convert substitutions/site distances to substitutions distances).
        :param reference_name: The name of the reference genome in the tree.
        :param reference_included: True if the reference genome is one of the leaves of the tree, False otherwise.
        :return: A new MutationTreeSamplesQuery object.
        """
        super().__init__(connection=connection, wrapped_query=wrapped_query, tree=tree)
        self._alignment_length = alignment_length
        self._reference_name = reference_name
        self._reference_included = reference_included

    @property
    def reference_name(self) -> str:
        """
        Gets the name of the reference genome.
        :return: The name of the reference genome.
        """
        return self._reference_name

    @property
    def reference_included(self) -> bool:
        """
        Whether or not the reference genome is included as a leaf of the tree.
        :return: True if the reference genome is a leaf of the tree, False otherwise.
        """
        return self._reference_included

    def _wrap_create(self, wrapped_query: SamplesQuery, universe_set: SampleSet = None) -> WrappedSamplesQuery:
        return MutationTreeSamplesQuery(connection=self._query_connection,
                                        wrapped_query=wrapped_query,
                                        tree=self._tree,
                                        alignment_length=self._alignment_length,
                                        reference_name=self._reference_name,
                                        reference_included=self._reference_included)

    def _within_distance_internal(self, data: Union[str, List[str], SamplesQuery, SampleSet], distance: float,
                                  units: str) -> SamplesQuery:
        if units == 'substitutions':
            distance_multiplier = self._alignment_length
        elif units == 'substitutions/site':
            distance_multiplier = 1
        else:
            raise Exception(f'Invalid units=[{units}]. Must be one of {self.DISTANCE_UNITS}')

        sample_names, query_infix = self._get_sample_names_query_infix_from_data(data)

        if len(sample_names) == 0:
            found_samples_set = SampleSet.create_empty()
        else:
            sample_name_ids_self = self._get_sample_name_ids(include_unknowns=True)
            found_samples_set = set()
            for sample_name in sample_names:
                sample_leaves = self._tree.get_leaves_by_name(sample_name)
                if len(sample_leaves) != 1:
                    raise Exception(
                        f'Invalid number of matching leaves for sample [{data}], leaves {sample_leaves}')

                sample_node = sample_leaves[0]

                for leaf in self._tree.iter_leaves():
                    if leaf.name not in sample_name_ids_self:
                        continue
                    sample_distance_to_other_sample = sample_node.get_distance(leaf) * distance_multiplier

                    if sample_distance_to_other_sample <= distance:
                        found_samples_set.add(sample_name_ids_self[leaf.name])
        found_samples = SampleSet(found_samples_set)
        return self.intersect(found_samples, f'within({distance} {units} of {query_infix})')

    def _can_handle_distance_units(self, units: str) -> bool:
        return units in self.DISTANCE_UNITS

    def _distance_units(self) -> List[str]:
        return self.DISTANCE_UNITS + self._wrapped_query._distance_units()
