import abc
from typing import Dict, List, Any

import pandas as pd

from genomics_data_index.api.query.features.FeaturesComparator import FeaturesComparator
from genomics_data_index.configuration.connector.DataIndexConnection import DataIndexConnection
from genomics_data_index.storage.SampleSet import SampleSet
from genomics_data_index.storage.model.db import FeatureSamples


class FeaturesFromIndexComparator(FeaturesComparator, abc.ABC):

    def __init__(self, connection: DataIndexConnection):
        super().__init__(connection=connection)

    def _create_summary_df(self, present_features: Dict[str, FeatureSamples],
                           present_samples: SampleSet) -> pd.DataFrame:
        data = []
        total = len(present_samples)
        for feature_id in present_features:
            feature = present_features[feature_id]
            samples_in_feature = present_samples.intersection(feature.sample_ids)
            sample_count = len(samples_in_feature)
            if sample_count > 0:
                data.append(self._create_feature_sample_count_row(feature_id,
                                                                  feature=feature,
                                                                  sample_count=sample_count,
                                                                  total=total))
        summary_df = pd.DataFrame(data,
                                  columns=[self.index_name] + self.summary_columns)
        return summary_df.set_index(self.index_name)

    def _create_summary_comparison_df(self, selected_samples: SampleSet,
                                      sample_categories: List[SampleSet],
                                      present_features: Dict[str, FeatureSamples],
                                      category_names: List[str] = None,
                                      compare_kind: str = 'percent') -> pd.DataFrame:
        data = []
        total = len(selected_samples)
        for feature_id in present_features:
            feature = present_features[feature_id]
            samples_in_feature = selected_samples.intersection(feature.sample_ids)

            for sample_category in sample_categories:
                samples_in_feature_in_cateogory = samples_in_feature.intersection(sample_category)
                sample_count = len(samples_in_feature_in_cateogory)
                if sample_count > 0:
                    data.append(self._create_feature_sample_count_row(feature_id,
                                                                      feature=feature,
                                                                      sample_count=sample_count,
                                                                      total=total))
        summary_df = pd.DataFrame(data,
                                  columns=[self.index_name] + self.summary_columns)
        return summary_df.set_index(self.index_name)

    @abc.abstractmethod
    def _create_feature_sample_count_row(self, feature_id: str, feature: FeatureSamples,
                                         sample_count: int, total: int) -> List[Any]:
        pass
