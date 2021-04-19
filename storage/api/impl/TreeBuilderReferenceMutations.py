from typing import Union, Iterable

from ete3 import Tree
import logging

from storage.api.TreeBuilder import TreeBuilder
from storage.connector.DataIndexConnection import DataIndexConnection
from storage.variant.SampleSet import SampleSet

logger = logging.getLogger(__name__)


class TreeBuilderReferenceMutations(TreeBuilder):
    TREE_METHODS = ['iqtree']

    def __init__(self, database_connection: DataIndexConnection, reference_name: str):
        super().__init__()
        self._database_connection = database_connection
        self._reference_name = reference_name

        if not database_connection.reference_service.exists_reference_genome(reference_name):
            raise Exception(f'Reference genome [{reference_name}] does not exist.')

    def build(self, samples_set: Union[SampleSet, Iterable[str]], method: str, **kwargs):
        if method == 'iqtree':
            return self._build_iqtree(samples_set=samples_set, **kwargs)
        else:
            raise Exception(f'Invalid method={method}, only support: {self.TREE_METHODS}')

    def _build_iqtree(self, samples_set: Union[SampleSet, Iterable[str]],
                      maxcores: int = 1,
                      align_type: str = 'full',
                      include_reference: bool = True,
                      extra_params: str = None) -> Tree:
        sample_service = self._database_connection.sample_service
        alignment_service = self._database_connection.alignment_service
        tree_service = self._database_connection.tree_service

        if isinstance(samples_set, SampleSet):
            samples_names = {s.name for s in sample_service.find_samples_by_ids(samples_set)}
        else:
            samples_names = set(samples_set)

        reference_samples = {s.name for s in sample_service.get_samples_associated_with_reference(self._reference_name)}

        tree_samples = samples_names - reference_samples

        if len(tree_samples) < len(samples_names):
            logger.warning(f'Building tree with only ({len(tree_samples)}/{len(reference_samples)} samples '
                           f'that are found on reference [{self._reference_name}]')

        alignment_data = alignment_service.construct_alignment(reference_name=self._reference_name,
                                                               samples=samples_names,
                                                               align_type=align_type,
                                                               include_reference=include_reference)
        tree_data, out = tree_service.build_tree(alignment_data, tree_build_type='iqtree',
                                                 num_cores=maxcores, align_type=align_type,
                                                 extra_params=extra_params)

        return tree_data
