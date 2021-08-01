import logging
from typing import Union, Iterable, Tuple

from ete3 import Tree

from genomics_data_index.api.query.TreeBuilder import TreeBuilder
from genomics_data_index.configuration.connector import DataIndexConnection
from genomics_data_index.storage.SampleSet import SampleSet

logger = logging.getLogger(__name__)


class TreeBuilderReferenceMutations(TreeBuilder):
    TREE_METHODS = ['iqtree']

    def __init__(self, database_connection: DataIndexConnection, reference_name: str):
        super().__init__()
        self._database_connection = database_connection
        self._reference_name = reference_name

        if not database_connection.reference_service.exists_reference_genome(reference_name):
            raise Exception(f'Reference genome [{reference_name}] does not exist.')

    def build(self, samples_set: Union[SampleSet, Iterable[str]],
              method: str = 'iqtree', include_reference=True, **kwargs) -> Tuple[Tree, int, SampleSet]:
        if method == 'iqtree':
            return self._build_iqtree(samples_set=samples_set, include_reference=include_reference,
                                      **kwargs)
        else:
            raise Exception(f'Invalid method={method}, only support: {self.TREE_METHODS}')

    def _build_iqtree(self, samples_set: Union[SampleSet, Iterable[str]],
                      ncores: int = 1,
                      align_type: str = 'full',
                      include_reference: bool = True,
                      extra_params: str = None) -> Tuple[Tree, int, SampleSet]:
        sample_service = self._database_connection.sample_service
        alignment_service = self._database_connection.alignment_service
        tree_service = self._database_connection.tree_service

        if isinstance(samples_set, SampleSet):
            samples_names = {s.name for s in sample_service.find_samples_by_ids(samples_set)}
        else:
            samples_names = set(samples_set)

        reference_samples = {s.name for s in sample_service.get_samples_associated_with_reference(self._reference_name)}

        tree_samples = samples_names.intersection(reference_samples)
        tree_samples_set = SampleSet(sample_service.find_sample_name_ids(tree_samples).values())

        if len(tree_samples) < len(samples_names):
            logger.warning(f'Building build_tree with only ({len(tree_samples)}/{len(reference_samples)}) samples '
                           f'that are found on reference [{self._reference_name}]')

        logger.info(f'Building build_tree using "iqtree" for {len(tree_samples_set)} samples')
        alignment_data = alignment_service.construct_alignment(reference_name=self._reference_name,
                                                               samples=samples_names,
                                                               align_type=align_type,
                                                               include_reference=include_reference)
        tree_data, out = tree_service.build_tree(alignment_data, tree_build_type='iqtree',
                                                 num_cores=ncores, align_type=align_type,
                                                 extra_params=extra_params)

        return tree_data, alignment_data.get_alignment_length(), tree_samples_set
