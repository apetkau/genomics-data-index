import pandas as pd

from genomics_data_index.configuration.connector.DataIndexConnection import DataIndexConnection
from genomics_data_index.storage.SampleSet import SampleSet
from genomics_data_index.api.query.features.FeaturesSummarizer import FeaturesSummarizer


class MutationFeaturesSummarizer(FeaturesSummarizer):

    def __init__(self, connection: DataIndexConnection, ignore_annotations: bool = False):
        super().__init__(connection=connection)
        self._ignore_annotations = ignore_annotations

    def summary(self, present_samples: SampleSet,
                unknown_samples: SampleSet,
                absent_samples: SampleSet,
                selection: str) -> pd.DataFrame:
        return self._summary_features_mutations(present_samples=present_samples,
                                                unknown_samples=unknown_samples,
                                                absent_samples=absent_samples,
                                                selection=selection)

    def summary_unique(self, samples_set: SampleSet, complement_set: SampleSet) -> pd.DataFrame:
        pass

    def _summary_features_mutations(self, present_samples: SampleSet,
                                    unknown_samples: SampleSet,
                                    absent_samples: SampleSet,
                                    selection: str = 'all',
                                    ncores: int = 1,
                                    batch_size: int = 500,
                                    mutation_type: str = 'all') -> pd.DataFrame:
        if selection not in self.FEATURES_SELECTIONS:
            raise Exception(f'selection=[{selection}] is unknown. Must be one of {self.FEATURES_SELECTIONS}')

        vs = self._connection.variation_service
        features_all_df = vs.count_mutations_in_sample_ids_dataframe(sample_ids=present_samples,
                                                                     ncores=ncores,
                                                                     batch_size=batch_size,
                                                                     mutation_type=mutation_type
                                                                     )
        features_all_df['Total'] = len(present_samples)
        features_all_df['Percent'] = 100 * (features_all_df['Count'] / features_all_df['Total'])

        if selection == 'all':
            features_results_df = features_all_df
        elif selection == 'unique':
            features_complement_df = self._summary_features_mutations(present_samples=absent_samples,
                                                  unknown_samples=unknown_samples,
                                                  absent_samples=present_samples,
                                                  selection='all')
            features_merged_df = features_all_df.merge(features_complement_df, left_index=True, right_index=True,
                                                       how='left', indicator=True, suffixes=('_x', '_y'))
            features_merged_df = features_merged_df[features_merged_df['_merge'] == 'left_only'].rename({
                'Sequence_x': 'Sequence',
                'Position_x': 'Position',
                'Deletion_x': 'Deletion',
                'Insertion_x': 'Insertion',
                'Count_x': 'Count',
                'Total_x': 'Total',
                'Percent_x': 'Percent',
            }, axis='columns')
            features_results_df = features_merged_df[['Sequence', 'Position', 'Deletion', 'Insertion',
                                                      'Count', 'Total', 'Percent']]
        else:
            raise Exception(f'selection=[{selection}] is unknown. Must be one of {self.FEATURES_SELECTIONS}')

        if not self._ignore_annotations:
            features_results_df = vs.append_mutation_annotations(features_results_df)

        return features_results_df
