from typing import List

import pandas as pd

from genomics_data_index.configuration.connector.DataIndexConnection import DataIndexConnection
from genomics_data_index.storage.SampleSet import SampleSet
from genomics_data_index.api.query.features.FeaturesSummarizer import FeaturesSummarizer


class MutationFeaturesSummarizer(FeaturesSummarizer):

    def __init__(self, connection: DataIndexConnection, ignore_annotations: bool = False,
                 ncores: int = 1, batch_size: int = 500, mutation_type: str = 'all'):
        super().__init__(connection=connection)
        self._ignore_annotations = ignore_annotations
        self._ncores = ncores
        self._batch_size = batch_size
        self._mutation_type = mutation_type

    def summary(self, sample_set: SampleSet) -> pd.DataFrame:
        vs = self._connection.variation_service
        features_df = vs.count_mutations_in_sample_ids_dataframe(sample_ids=sample_set,
                                                                 ncores=self._ncores,
                                                                 batch_size=self._batch_size,
                                                                 mutation_type=self._mutation_type
                                                                 )
        features_df['Total'] = len(sample_set)
        features_df['Percent'] = 100 * (features_df['Count'] / features_df['Total'])

        return self._join_additional_columns(features_df)

    def _join_additional_columns(self, features_df: pd.DataFrame) -> pd.DataFrame:
        if not self._ignore_annotations:
            return self._connection.variation_service.append_mutation_annotations(features_df)
        else:
            return features_df

    @property
    def summary_columns(self) -> List[str]:
        return ['Sequence', 'Position', 'Deletion', 'Insertion', 'Count', 'Total', 'Percent']
