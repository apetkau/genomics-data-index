from pathlib import Path
from typing import List, Dict

import pandas as pd

from storage.variant.service.QueryService import QueryService
from storage.variant.service.ReferenceService import ReferenceService
from storage.variant.service.TreeService import TreeService
from storage.variant.service.SampleService import SampleService

from storage.variant.service.QueryService import QueryFeature


class QueryFeatureMutation(QueryFeature):

    def __init__(self, spdi: str):
        super().__init__()
        self._spdi = spdi

    @property
    def spdi(self):
        return self._spdi


class MutationQueryService(QueryService):

    def __init__(self, tree_service: TreeService, reference_service: ReferenceService,
                 sample_service: SampleService):
        super().__init__()
        self._tree_service = tree_service
        self._reference_service = reference_service
        self._sample_service = sample_service

    def _find_matches_internal(self, sample_names: List[str], distance_threshold: float = None) -> pd.DataFrame:
        sample_distances = []
        for sample_name in sample_names:
            reference_genomes = self._reference_service.find_references_for_sample(sample_name)

            for reference_genome in reference_genomes:
                tree = reference_genome.tree

                sample_leaves = tree.get_leaves_by_name(sample_name)
                if len(sample_leaves) != 1:
                    raise Exception(
                        f'Invalid number of matching leaves for sample [{sample_name}], leaves {sample_leaves}')

                sample_node = sample_leaves[0]

                leaves = tree.get_leaves()
                for leaf in leaves:
                    if leaf.name == sample_node.name:
                        continue
                    distance = sample_node.get_distance(leaf)
                    align_length = reference_genome.tree_alignment_length
                    sample_distances.append([reference_genome.name, sample_name, leaf.name,
                                             f'{distance * align_length:0.2f}', distance, align_length])

        matches_df = pd.DataFrame(data=sample_distances, columns=[
            'Reference Genome',
            'Sample A',
            'Sample B',
            'Distance',
            'Distance (subs/site)',
            'SNV Alignment Length',
        ])
        matches_df['Distance'] = pd.to_numeric(matches_df['Distance'])
        matches_df['Distance (subs/site)'] = pd.to_numeric(matches_df['Distance (subs/site)'])
        matches_df['SNV Alignment Length'] = pd.to_numeric(matches_df['SNV Alignment Length'])
        matches_df = matches_df.sort_values(['Sample A', 'Distance (subs/site)'], ascending=True)

        if distance_threshold is not None:
            matches_df = matches_df.loc[:, matches_df['Distance (subs/site)'] <= distance_threshold]

        return matches_df

    def _find_by_features_internal(self, features: List[QueryFeature]) -> pd.DataFrame:
        for feature in features:
            if not isinstance(feature, QueryFeatureMutation):
                raise Exception(f'feature=[{feature}] is not of type QueryFeatureMutation')

        variation_ids = [f.spdi for f in features]
        variation_samples = self._sample_service.find_samples_by_variation_ids(variation_ids)

        data = []
        for vid in variation_samples:
            for sample in variation_samples[vid]:
                data.append([vid, sample.name, sample.id])

        return pd.DataFrame(data=data, columns=['Feature', 'Sample Name', 'Sample ID'])

    def _find_matches_genome_files_internal(self, sample_reads: Dict[str, List[Path]],
                                            distance_threshold: float = None) -> pd.DataFrame:
        raise Exception('Method not implemented')

    def _pairwise_distance_internal(self, samples: List[str]) -> pd.DataFrame:
        raise Exception('Method not implemented')

    def _differences_between_genomes_internal(self, sample1: str, sample2: str):
        raise Exception('Method not implemented')

    def get_data_type(self) -> str:
        return 'SNV'
