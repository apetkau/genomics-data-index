from pathlib import Path
from typing import List, Dict

import pandas as pd
from pandas.api.types import CategoricalDtype

from storage.variant.service.QueryService import QueryFeature
from storage.variant.service.QueryService import QueryService
from storage.variant.service.ReferenceService import ReferenceService
from storage.variant.service.SampleService import SampleService
from storage.variant.service.TreeService import TreeService


class QueryFeatureMutation(QueryFeature):

    def __init__(self, spdi: str):
        super().__init__()
        self._spdi = spdi

        seq, pos, ref, alt = self._spdi.split(':')
        self._seq_name = seq
        self._pos = int(pos)
        self._ref = ref
        self._alt = alt

    @property
    def spdi(self):
        return self._spdi

    @property
    def sequence_name(self):
        return self._seq_name

    @property
    def position(self):
        return self._pos

    @property
    def start(self):
        return self.position

    @property
    def stop(self):
        return self.position + len(self.ref)

    @property
    def ref(self):
        return self._ref

    @property
    def alt(self):
        return self._alt


class MutationQueryService(QueryService):

    def __init__(self, reference_service: ReferenceService,
                 sample_service: SampleService,
                 tree_service: TreeService):
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

    def _find_by_features_internal(self, features: List[QueryFeature], include_unknown: bool) -> pd.DataFrame:
        for feature in features:
            if not isinstance(feature, QueryFeatureMutation):
                raise Exception(f'feature=[{feature}] is not of type QueryFeatureMutation')

        variation_ids = [f.spdi for f in features]
        variation_samples = self._sample_service.find_samples_by_variation_ids(variation_ids)

        data = []
        for vid in variation_samples:
            for sample in variation_samples[vid]:
                data.append([vid, sample.name, sample.id, 'Present'])

        if include_unknown:
            for feature in features:
                reference = self._reference_service.find_reference_for_sequence(feature.sequence_name)

                # TODO: I'm doing linear search over all samples, this could be improved
                for sample_variation in reference.sample_nucleotide_variation:
                    masked_regions = sample_variation.masked_regions
                    if masked_regions.overlaps_range(feature.sequence_name, feature.start, feature.stop):
                        data.append([feature.spdi, sample_variation.sample.name, sample_variation.sample.id, 'Unknown'])

        return pd.DataFrame(data=data, columns=[
            'Feature', 'Sample Name', 'Sample ID', 'Status']).sort_values(['Feature', 'Sample Name'])

    def _find_matches_genome_files_internal(self, sample_reads: Dict[str, List[Path]],
                                            distance_threshold: float = None) -> pd.DataFrame:
        raise Exception('Method not implemented')

    def _pairwise_distance_internal(self, samples: List[str]) -> pd.DataFrame:
        raise Exception('Method not implemented')

    def _differences_between_genomes_internal(self, sample1: str, sample2: str):
        raise Exception('Method not implemented')

    def get_data_type(self) -> str:
        return 'mutation'


class MutationQuerySummaries:

    def __init__(self):
        pass

    def find_by_features_summary(self, find_by_features_results: pd.DataFrame,
                                 sample_counts: Dict[str, int]) -> pd.DataFrame:
        summary_df = find_by_features_results.astype({
            'Status': CategoricalDtype(categories=['Present', 'Absent', 'Unknown'])
        })

        summary_df = summary_df.groupby(['Feature', 'Status']).agg({'Sample Name': 'count'})
        summary_df = summary_df.unstack().fillna(0).astype('int64')
        summary_df.columns = list(summary_df.columns.droplevel())

        summary_df['Total'] = summary_df.apply(lambda x: sample_counts[x.name], axis='columns')
        summary_df['Absent'] = summary_df['Total'] - (summary_df['Present'] + summary_df['Unknown'])
        summary_df['% Present'] = 100 * summary_df['Present'] / summary_df['Total']
        summary_df['% Absent'] = 100 * summary_df['Absent'] / summary_df['Total']
        summary_df['% Unknown'] = 100 * summary_df['Unknown'] / summary_df['Total']
        summary_df.columns.name = None

        return summary_df
