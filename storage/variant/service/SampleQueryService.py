import pandas as pd

from storage.variant.service.TreeService import TreeService
from storage.variant.service.ReferenceService import ReferenceService


class SampleQueryService:

    def __init__(self, tree_service: TreeService, reference_service: ReferenceService):
        self._tree_service = tree_service
        self._reference_service = reference_service

    def find_matches(self, sample_name: str) -> pd.DataFrame:
        reference_genomes = self._reference_service.find_references_for_sample(sample_name)

        sample_distances = []
        for reference_genome in reference_genomes:
            tree = reference_genome.tree

            sample_leaves = tree.get_leaves_by_name(sample_name)
            if len(sample_leaves) != 1:
                raise Exception(f'Invalid number of matching leaves for sample [{sample_name}], leaves {sample_leaves}')

            sample_node = sample_leaves[0]

            leaves = tree.get_leaves()
            for leaf in leaves:
                if leaf.name == sample_node.name:
                    continue
                distance = sample_node.get_distance(leaf)
                align_length = reference_genome.tree_alignment_length
                sample_distances.append([reference_genome.name, sample_name, leaf.name,
                                         distance, f'{distance * align_length:0.2f}', align_length])

        matches_df = pd.DataFrame(data=sample_distances, columns=[
            'Reference Genome',
            'Sample A',
            'Sample B',
            'Distance (subs/site)',
            'Distance (subs)',
            'Alignment Length',
        ])
        matches_df['Distance (subs/site)'] = pd.to_numeric(matches_df['Distance (subs/site)'])
        matches_df['Distance (subs)'] = pd.to_numeric(matches_df['Distance (subs)'])
        matches_df['Alignment Length'] = pd.to_numeric(matches_df['Alignment Length'])
        return matches_df.sort_values('Distance (subs/site)', ascending=True)