from typing import List

import pandas as pd

from genomics_data_index.api.query.features.FeaturesComparator import FeaturesComparator
from genomics_data_index.configuration.connector.DataIndexConnection import DataIndexConnection
from genomics_data_index.storage.SampleSet import SampleSet


class MutationFeaturesComparator(FeaturesComparator):

    def __init__(self, connection: DataIndexConnection, ignore_annotations: bool = False,
                 include_unknown_samples: bool = True,
                 ncores: int = 1, batch_size: int = 500, mutation_type: str = 'all'):
        super().__init__(connection=connection, include_unknown_samples=include_unknown_samples)
        self._ignore_annotations = ignore_annotations
        self._ncores = ncores
        self._batch_size = batch_size
        self._mutation_type = mutation_type

    @property
    def summary_columns(self) -> List[str]:
        return self.feature_id_columns + ['Count', 'Total', 'Percent']

    @property
    def feature_id_columns(self) -> List[str]:
        return ['Sequence', 'Position', 'Deletion', 'Insertion']

    @property
    def index_name(self) -> str:
        return 'Mutation'

    def features_comparison(self, selected_samples: SampleSet,
                            sample_categories: List[SampleSet],
                            category_prefixes: List[str] = None,
                            category_samples_threshold: int = None,
                            unit: str = 'percent') -> pd.DataFrame:
        raise NotImplementedError()

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
