import abc
from typing import Dict, List, Any

import pandas as pd

from genomics_data_index.api.query.features.FeaturesComparator import FeaturesComparator
from genomics_data_index.configuration.connector.DataIndexConnection import DataIndexConnection
from genomics_data_index.storage.SampleSet import SampleSet
from genomics_data_index.storage.model.db import FeatureSamples


class FeatureSamplesSummarizer(abc.ABC):

    def __init__(self):
        pass

    @abc.abstractmethod
    def summary_data(self, samples: SampleSet, total: int) -> List[Any]:
        pass

    @abc.abstractmethod
    def summary_names(self) -> List[str]:
        pass


class FeatureSamplesSingleCategorySummarizer(FeatureSamplesSummarizer):

    SUMMARY_NAMES = ['Count', 'Total', 'Percent']

    def __init__(self):
        super().__init__()

    def summary_data(self, samples: SampleSet, total: int) -> List[Any]:
        sample_count = len(samples)
        percent = (sample_count / total) * 100
        return [sample_count, total, percent]

    def summary_names(self) -> List[str]:
        return self.SUMMARY_NAMES


class FeaturesFromIndexComparator(FeaturesComparator, abc.ABC):

    def __init__(self, connection: DataIndexConnection):
        super().__init__(connection=connection)

    def _create_summary_df(self, present_features: Dict[str, FeatureSamples],
                           present_samples: SampleSet,
                           feature_samples_summarizer: FeatureSamplesSummarizer) -> pd.DataFrame:
        data = []
        total = len(present_samples)
        for feature_id in present_features:
            feature = present_features[feature_id]
            samples_in_feature = present_samples.intersection(feature.sample_ids)
            sample_count = len(samples_in_feature)
            if sample_count > 0:
                data.append(self._create_feature_sample_count_row(feature_id,
                                                                  feature=feature,
                                                                  feature_samples=samples_in_feature,
                                                                  total=total,
                                                                  feature_samples_summarizer=feature_samples_summarizer))
        summary_df = pd.DataFrame(data,
                                  columns=[self.index_name] + self.feature_id_columns + feature_samples_summarizer.summary_names())
        return summary_df.set_index(self.index_name)

    # def _create_summary_comparison_df(self, selected_samples: SampleSet,
    #                                   sample_categories: List[SampleSet],
    #                                   present_features: Dict[str, FeatureSamples],
    #                                   category_names: List[str] = None,
    #                                   compare_kind: str = 'percent') -> pd.DataFrame:
    #     data = []
    #     total = len(selected_samples)
    #     for feature_id in present_features:
    #         feature = present_features[feature_id]
    #         samples_in_feature = selected_samples.intersection(feature.sample_ids)
    #
    #         for sample_category in sample_categories:
    #             samples_in_feature_in_cateogory = samples_in_feature.intersection(sample_category)
    #             sample_count = len(samples_in_feature_in_cateogory)
    #             if sample_count > 0:
    #                 data.append(self._create_feature_sample_count_row(feature_id,
    #                                                                   feature=feature,
    #                                                                   sample_count=sample_count,
    #                                                                   total=total))
    #     summary_df = pd.DataFrame(data,
    #                               columns=[self.index_name] + self.summary_columns)
    #     return summary_df.set_index(self.index_name)

    @abc.abstractmethod
    def _create_feature_sample_count_row(self, feature_id: str, feature: FeatureSamples,
                                         feature_samples: SampleSet,
                                         total: int,
                                         feature_samples_summarizer: FeatureSamplesSummarizer) -> List[Any]:
        pass
