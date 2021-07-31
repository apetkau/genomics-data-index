from typing import List, Any, cast

import pandas as pd

from genomics_data_index.api.query.features.FeaturesFromIndexSummarizer import FeaturesFromIndexSummarizer
from genomics_data_index.configuration.connector.DataIndexConnection import DataIndexConnection
from genomics_data_index.storage.SampleSet import SampleSet
from genomics_data_index.storage.model.QueryFeatureMutationSPDI import QueryFeatureMutationSPDI
from genomics_data_index.storage.model.db import FeatureSamples
from genomics_data_index.storage.model.db import NucleotideVariantsSamples
from genomics_data_index.storage.service.VariationService import VariationService


class MutationFeaturesFromIndexSummarizer(FeaturesFromIndexSummarizer):

    def __init__(self, connection: DataIndexConnection, ignore_annotations: bool = False,
                 include_present: bool = True, include_unknown: bool = False,
                 mutation_type: str = 'all', id_type: str = 'spdi_ref'):
        super().__init__(connection=connection)
        self._ignore_annotations = ignore_annotations
        self._mutation_type = mutation_type
        self._include_present = include_present
        self._include_unknown = include_unknown
        self._id_type = id_type

    @property
    def summary_columns(self) -> List[str]:
        return ['Sequence', 'Position', 'Deletion', 'Insertion', 'Count', 'Total', 'Percent']

    @property
    def index_name(self) -> str:
        return 'Mutation'

    def _create_feature_sample_count_row(self, feature_id: str, feature: FeatureSamples,
                                         sample_count: int, total: int) -> List[Any]:
        if isinstance(feature, NucleotideVariantsSamples):
            feature = cast(NucleotideVariantsSamples, feature)
            percent = (sample_count / total) * 100
            query_feature_id = QueryFeatureMutationSPDI(
                feature_id)  # Use this since db deletion is an int (not sequence)
            return [feature_id, feature.sequence, feature.position, query_feature_id.deletion, feature.insertion,
                    sample_count, total, percent]
        else:
            raise Exception(f'feature={feature} is not of type {NucleotideVariantsSamples.__name__}')

    def _join_additional_columns(self, features_df: pd.DataFrame) -> pd.DataFrame:
        if not self._ignore_annotations:
            return self._connection.variation_service.append_mutation_annotations(features_df)
        else:
            return features_df

    def summary(self, sample_set: SampleSet) -> pd.DataFrame:
        variation_service: VariationService = self._connection.variation_service
        present_features = variation_service.get_features(include_present=self._include_present,
                                                          include_unknown=self._include_unknown,
                                                          id_type=self._id_type)

        features_df = self._create_summary_df(present_features, present_samples=sample_set)

        return self._join_additional_columns(features_df)
