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


class FeatureSamplesMultipleCategorySummarizer(FeatureSamplesSummarizer):

    def __init__(self, sample_categories: List[SampleSet], category_prefixes: List[str] = None,
                 compare_kind: str = 'percent'):
        super().__init__()
        self._sample_categories = sample_categories

        if category_prefixes is None:
            category_prefixes = [f'Category{x + 1}' for x in range(len(sample_categories))]
        elif not isinstance(category_prefixes, list):
            raise Exception(f'category_names={category_prefixes} must be a list or None')
        elif len(category_prefixes) != len(sample_categories):
            raise Exception(f'sample_categories has {len(sample_categories)} elements but, '
                            f'category_names has {len(category_prefixes)} elements. These must be the same size.')

        # I convert to string in case someone passes a list that's not composed of strings
        self._summary_names = ['Total'] + [str(c) for c in category_prefixes]

        if compare_kind not in ['percent', 'count']:
            raise Exception(f'compare_kind={compare_kind} must be one of "percent" or "count"')
        elif compare_kind == 'percent':
            self._use_percent = True
        else:
            self._use_percent = False

    def summary_data(self, samples: SampleSet, total: int) -> List[Any]:
        data = [total]
        for sample_category in self._sample_categories:
            samples_in_category = samples.intersection(sample_category)
            category_count = len(samples_in_category)
            if self._use_percent:
                category_percent = (category_count / total) * 100
                data.append(category_percent)
            else:
                data.append(category_count)
        return data

    def summary_names(self) -> List[str]:
        return self._summary_names


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

    def summary(self, sample_set: SampleSet) -> pd.DataFrame:
        samples_summarizer = FeatureSamplesSingleCategorySummarizer()
        return self._do_summary(sample_set=sample_set, feature_samples_summarizer=samples_summarizer)

    def features_comparison(self, selected_samples: SampleSet,
                            sample_categories: List[SampleSet],
                            category_prefixes: List[str] = None,
                            compare_kind: str = 'percent') -> pd.DataFrame:
        samples_summarizer = FeatureSamplesMultipleCategorySummarizer(sample_categories=sample_categories,
                                                                      category_prefixes=category_prefixes,
                                                                      compare_kind=compare_kind)
        return self._do_summary(sample_set=selected_samples, feature_samples_summarizer=samples_summarizer)

    @abc.abstractmethod
    def _do_summary(self, sample_set: SampleSet, feature_samples_summarizer: FeatureSamplesSummarizer) -> pd.DataFrame:
        pass

    @abc.abstractmethod
    def _create_feature_sample_count_row(self, feature_id: str, feature: FeatureSamples,
                                         feature_samples: SampleSet,
                                         total: int,
                                         feature_samples_summarizer: FeatureSamplesSummarizer) -> List[Any]:
        pass
