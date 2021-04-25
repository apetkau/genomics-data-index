from pathlib import Path
from typing import List, Dict, Any, Set

import pandas as pd

from genomics_data_index.variant.model.QueryFeature import QueryFeature
from genomics_data_index.variant.model.QueryFeatureMutation import QueryFeatureMutation
from genomics_data_index.variant.service.FullFeatureQueryService import FullFeatureQueryService
from genomics_data_index.variant.service.ReferenceService import ReferenceService
from genomics_data_index.variant.service.SampleService import SampleService
from genomics_data_index.variant.service.TreeService import TreeService


class MutationQueryService(FullFeatureQueryService):

    def __init__(self, reference_service: ReferenceService,
                 sample_service: SampleService,
                 tree_service: TreeService):
        super().__init__(sample_service)
        self._tree_service = tree_service
        self._reference_service = reference_service

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
            'Query',
            'Match',
            'Distance',
            'Distance (subs/site)',
            'SNV Alignment Length',
        ])
        matches_df['Distance'] = pd.to_numeric(matches_df['Distance'])
        matches_df['Distance (subs/site)'] = pd.to_numeric(matches_df['Distance (subs/site)'])
        matches_df['SNV Alignment Length'] = pd.to_numeric(matches_df['SNV Alignment Length'])
        matches_df = matches_df.sort_values(['Query', 'Distance (subs/site)'], ascending=True)

        if distance_threshold is not None:
            matches_df = matches_df.loc[:, matches_df['Distance (subs/site)'] <= distance_threshold]

        return matches_df

    def get_correct_query_feature(self) -> Any:
        return QueryFeatureMutation

    def expand_feature(self, feature: QueryFeature):
        raise Exception('Not implemented')

    def _get_feature_scope_sample_counts(self, feature_scopes: Set[str]) -> Dict[str, int]:
        sequence_reference_map = {n: self._reference_service.find_reference_for_sequence(n)
                                  for n in feature_scopes}
        sequence_sample_counts = {
            n: self._sample_service.count_samples_associated_with_reference(
                sequence_reference_map[n].name) for n in sequence_reference_map
        }

        return sequence_sample_counts

    def _get_unknown_features(self, features: List[QueryFeature]) -> pd.DataFrame:
        data = []
        for feature in features:
            reference = self._reference_service.find_reference_for_sequence(feature.scope)

            # TODO: I'm doing linear search over all samples, this could be improved
            for sample_variation in reference.sample_nucleotide_variation:
                masked_regions = sample_variation.masked_regions
                if masked_regions.overlaps_range(feature.scope, feature.start, feature.stop):
                    data.append([feature.id, sample_variation.sample.name, sample_variation.sample.id, 'Unknown'])

        return pd.DataFrame(data, columns=['Feature', 'Sample Name', 'Sample ID', 'Status'])

    def _find_matches_genome_files_internal(self, sample_reads: Dict[str, List[Path]],
                                            distance_threshold: float = None) -> pd.DataFrame:
        raise Exception('Method not implemented')

    def _pairwise_distance_internal(self, samples: List[str]) -> pd.DataFrame:
        raise Exception('Method not implemented')

    def _differences_between_genomes_internal(self, sample1: str, sample2: str):
        raise Exception('Method not implemented')

    def get_data_type(self) -> str:
        return 'mutation'
