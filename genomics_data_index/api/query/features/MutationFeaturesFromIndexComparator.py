from typing import List, Any, cast, Optional

import pandas as pd

from genomics_data_index.api.query.features.FeaturesFromIndexComparator import FeatureSamplesSingleCategorySummarizer
from genomics_data_index.api.query.features.FeaturesFromIndexComparator import FeatureSamplesSummarizer
from genomics_data_index.api.query.features.FeaturesFromIndexComparator import FeaturesFromIndexComparator
from genomics_data_index.configuration.connector.DataIndexConnection import DataIndexConnection
from genomics_data_index.storage.SampleSet import SampleSet
from genomics_data_index.storage.model.QueryFeatureMutationSPDI import QueryFeatureMutationSPDI
from genomics_data_index.storage.model.db import FeatureSamples
from genomics_data_index.storage.model.db import NucleotideVariantsSamples
from genomics_data_index.storage.service.VariationService import VariationService


class MutationFeaturesFromIndexComparator(FeaturesFromIndexComparator):

    def __init__(self, connection: DataIndexConnection, ignore_annotations: bool = False,
                 include_present: bool = True, include_unknown: bool = False,
                 include_unknown_samples: bool = True,
                 mutation_type: str = 'all', id_type: str = 'spdi_ref'):
        super().__init__(connection=connection, include_unknown_samples=include_unknown_samples)
        self._ignore_annotations = ignore_annotations
        self._mutation_type = mutation_type
        self._include_present = include_present
        self._include_unknown = include_unknown
        self._id_type = id_type

    @property
    def summary_columns(self) -> List[str]:
        return self.feature_id_columns + FeatureSamplesSingleCategorySummarizer.SUMMARY_NAMES

    @property
    def feature_id_columns(self) -> List[str]:
        return ['Sequence', 'Position', 'Deletion', 'Insertion', 'Type']

    @property
    def index_name(self) -> str:
        return 'Mutation'

    def _create_feature_sample_count_row(self, feature_id: str, feature: FeatureSamples,
                                         feature_samples: SampleSet,
                                         feature_samples_unknown: Optional[SampleSet],
                                         total: int,
                                         feature_samples_summarizer: FeatureSamplesSummarizer) -> List[Any]:
        if isinstance(feature, NucleotideVariantsSamples):
            feature = cast(NucleotideVariantsSamples, feature)
            summary_data = feature_samples_summarizer.summary_data(samples=feature_samples,
                                                                   unknown_samples=feature_samples_unknown,
                                                                   total=total)
            query_feature_id = QueryFeatureMutationSPDI(
                feature_id)  # Use this since db deletion is an int (not sequence)
            return [feature_id, feature.sequence, feature.position,
                    query_feature_id.deletion, feature.insertion, feature.var_type] + summary_data
        else:
            raise Exception(f'feature={feature} is not of type {NucleotideVariantsSamples.__name__}')

    def _join_additional_columns(self, features_df: pd.DataFrame) -> pd.DataFrame:
        if not self._ignore_annotations:
            return self._connection.variation_service.append_mutation_annotations(features_df)
        else:
            return features_df

    def _do_summary(self, sample_set: SampleSet, feature_samples_summarizer: FeatureSamplesSummarizer) -> pd.DataFrame:
        variation_service: VariationService = self._connection.variation_service
        present_features = variation_service.get_features(include_present=self._include_present,
                                                          include_unknown=self._include_unknown,
                                                          id_type=self._id_type)

        features_df = self._create_summary_df(present_features, present_samples=sample_set,
                                              feature_samples_summarizer=feature_samples_summarizer)

        return self._join_additional_columns(features_df)
