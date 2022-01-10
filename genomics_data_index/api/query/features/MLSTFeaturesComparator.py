from typing import List, Any, Optional

import pandas as pd

from genomics_data_index.api.query.features.FeaturesFromIndexComparator import FeatureSamplesSingleCategorySummarizer
from genomics_data_index.api.query.features.FeaturesFromIndexComparator import FeatureSamplesSummarizer
from genomics_data_index.api.query.features.FeaturesFromIndexComparator import FeaturesFromIndexComparator
from genomics_data_index.configuration.connector.DataIndexConnection import DataIndexConnection
from genomics_data_index.storage.SampleSet import SampleSet
from genomics_data_index.storage.model.db import FeatureSamples, MLSTAllelesSamples
from genomics_data_index.storage.service.MLSTService import MLSTService


class MLSTFeaturesComparator(FeaturesFromIndexComparator):

    def __init__(self, connection: DataIndexConnection,
                 scheme: str = None,
                 locus: str = None,
                 include_unknown_samples: bool = True,
                 include_unknown_no_present_samples=False,
                 include_present: bool = True,
                 include_unknown: bool = False):
        super().__init__(connection=connection, include_unknown_samples=include_unknown_samples,
                         include_unknown_no_present_samples=include_unknown_no_present_samples)
        self._scheme = scheme
        self._locus = locus
        self._include_present = include_present
        self._include_unknown = include_unknown

    @property
    def summary_columns(self) -> List[str]:
        return self.feature_id_columns + FeatureSamplesSingleCategorySummarizer.SUMMARY_NAMES

    @property
    def feature_id_columns(self) -> List[str]:
        return ['Scheme', 'Locus', 'Allele']

    @property
    def index_name(self) -> str:
        return 'MLST Feature'

    def _create_feature_sample_count_row(self, feature_id: str, feature: FeatureSamples,
                                         feature_samples: SampleSet,
                                         feature_samples_unknown: Optional[SampleSet],
                                         total: int,
                                         feature_samples_summarizer: FeatureSamplesSummarizer) -> List[Any]:
        if isinstance(feature, MLSTAllelesSamples):
            summary_data = feature_samples_summarizer.summary_data(samples=feature_samples, total=total,
                                                                   unknown_samples=feature_samples_unknown)
            return [feature.query_id, feature.scheme, feature.locus, feature.allele] + summary_data
        else:
            raise Exception(f'feature={feature} is not of type {MLSTAllelesSamples.__name__}')

    def _do_summary(self, sample_set: SampleSet, feature_samples_summarizer: FeatureSamplesSummarizer) -> pd.DataFrame:
        mlst_service: MLSTService = self._connection.mlst_service
        mlst_present_features = mlst_service.get_features(scheme=self._scheme,
                                                          locus=self._locus,
                                                          include_present=self._include_present,
                                                          include_unknown=self._include_unknown)
        return self._create_summary_df(mlst_present_features, present_samples=sample_set,
                                       feature_samples_summarizer=feature_samples_summarizer)
