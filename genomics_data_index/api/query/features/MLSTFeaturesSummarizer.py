from typing import List, cast, Any

import pandas as pd

from genomics_data_index.api.query.features.FeaturesFromIndexSummarizer import FeaturesFromIndexSummarizer
from genomics_data_index.configuration.connector.DataIndexConnection import DataIndexConnection
from genomics_data_index.storage.service.MLSTService import MLSTService
from genomics_data_index.storage.SampleSet import SampleSet
from genomics_data_index.storage.model.db import FeatureSamples, MLSTAllelesSamples


class MLSTFeaturesSummarizer(FeaturesFromIndexSummarizer):

    def __init__(self, connection: DataIndexConnection,
                 scheme: str = None,
                 locus: str = None,
                 include_present: bool = True,
                 include_unknown: bool = False):
        super().__init__(connection=connection)
        self._scheme = scheme
        self._locus = locus
        self._include_present = include_present
        self._include_unknown = include_unknown

    @property
    def summary_columns(self) -> List[str]:
        return ['Scheme', 'Locus', 'Allele', 'Count', 'Total', 'Percent']

    def _create_feature_sample_count_row(self, feature: FeatureSamples,
                                         sample_count: int, total: int) -> List[Any]:
        if isinstance(feature, MLSTAllelesSamples):
            percent = (sample_count / total) * 100
            return [feature.query_id, feature.scheme, feature.locus, feature.allele,
                    sample_count, total, percent]
        else:
            raise Exception(f'feature={feature} is not of type {MLSTAllelesSamples.__name__}')

    def summary(self, sample_set: SampleSet) -> pd.DataFrame:
        mlst_service: MLSTService = self._connection.mlst_service

        mlst_present_features = mlst_service.get_features(scheme=self._scheme,
                                                          locus=self._locus,
                                                          include_present=self._include_present,
                                                          include_unknown=self._include_unknown)
        data = self._create_feature_sample_count_rows(mlst_present_features, present_samples=sample_set)
        summary_df = pd.DataFrame(data,
                                  columns=['MLST Feature', 'Scheme', 'Locus', 'Allele', 'Count',
                                           'Total', 'Percent'])
        return summary_df.set_index('MLST Feature')
