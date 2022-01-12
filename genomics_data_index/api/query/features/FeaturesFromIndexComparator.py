import abc
from typing import Dict, List, Any, Optional

import logging
import pandas as pd

from genomics_data_index.api.query.features.FeaturesComparator import FeaturesComparator
from genomics_data_index.configuration.connector.DataIndexConnection import DataIndexConnection
from genomics_data_index.storage.SampleSet import SampleSet
from genomics_data_index.storage.model.QueryFeature import QueryFeature
from genomics_data_index.storage.model.QueryFeatureFactory import QueryFeatureFactory
from genomics_data_index.storage.model.db import FeatureSamples
from genomics_data_index.storage.service.SampleService import SampleService


logger = logging.getLogger(__name__)


class FeatureSamplesSummarizer(abc.ABC):

    def __init__(self):
        pass

    @abc.abstractmethod
    def summary_data(self, samples: SampleSet, unknown_samples: Optional[SampleSet], total: int) -> List[Any]:
        pass

    @abc.abstractmethod
    def summary_names(self) -> List[str]:
        pass


class FeatureSamplesSingleCategorySummarizer(FeatureSamplesSummarizer):
    SUMMARY_NAMES = ['Count', 'Unknown Count', 'Present and Unknown Count',
                     'Total', 'Percent', 'Unknown Percent', 'Present and Unknown Percent']

    def __init__(self):
        super().__init__()

    def summary_data(self, samples: SampleSet, unknown_samples: Optional[SampleSet], total: int) -> List[Any]:
        sample_count = len(samples)

        if unknown_samples is not None:
            unknown_sample_count = len(unknown_samples)
            present_and_unknown_sample_count = sample_count + unknown_sample_count
            unknown_percent = (unknown_sample_count / total) * 100
            present_and_unknown_percent = (present_and_unknown_sample_count / total) * 100
        else:
            unknown_sample_count = pd.NA
            present_and_unknown_sample_count = pd.NA
            unknown_percent = pd.NA
            present_and_unknown_percent = pd.NA

        percent = (sample_count / total) * 100
        return [sample_count, unknown_sample_count, present_and_unknown_sample_count, total,
                percent, unknown_percent, present_and_unknown_percent]

    def summary_names(self) -> List[str]:
        return self.SUMMARY_NAMES


class FeatureSamplesMultipleCategorySummarizer(FeatureSamplesSummarizer):

    def __init__(self, sample_categories: List[SampleSet], category_prefixes: List[str] = None,
                 compare_kind: str = 'percent'):
        super().__init__()
        sample_categories_totals = [len(c) for c in sample_categories]
        self._sample_categories_and_totals = list(zip(sample_categories, sample_categories_totals))

        self._use_count = False
        self._factor = None
        if compare_kind not in ['proportion', 'percent', 'count']:
            raise Exception(f'compare_kind={compare_kind} must be one of "percent", "proportion", or "count"')
        elif compare_kind == 'proportion':
            self._factor = 1
        elif compare_kind == 'percent':
            self._factor = 100
        else:
            self._use_count = True

        if category_prefixes is None:
            category_prefixes = [f'Category{x + 1}' for x in range(len(sample_categories))]
        elif not isinstance(category_prefixes, list):
            raise Exception(f'category_names={category_prefixes} must be a list or None')
        elif len(category_prefixes) != len(sample_categories):
            raise Exception(f'sample_categories has {len(sample_categories)} elements but, '
                            f'category_names has {len(category_prefixes)} elements. These must be the same size.')
        else:
            # I convert to string in case someone passes a list that's not composed of strings
            category_prefixes = [str(c) for c in category_prefixes]

        self._summary_names = ['Total'] \
                              + [f'{c}_{compare_kind}' for c in category_prefixes] \
                              + [f'{c}_total' for c in category_prefixes]

    def summary_data(self, samples: SampleSet, unknown_samples: Optional[SampleSet], total: int) -> List[Any]:
        data = [total]
        for sample_category, sample_category_total in self._sample_categories_and_totals:
            samples_in_category = samples.intersection(sample_category)
            category_count = len(samples_in_category)

            if unknown_samples is not None:
                unknown_samples_in_category = unknown_samples.intersection(sample_category)
                unknown_sample_count = len(unknown_samples_in_category)
                unknown_percent = (unknown_sample_count / total) * 100
            else:
                unknown_sample_count = pd.NA
                unknown_percent = pd.NA

            if self._use_count:
                data.append(category_count)
            else:
                if sample_category_total > 0:
                    category_percent = (category_count / sample_category_total) * self._factor
                else:
                    category_percent = pd.NA
                data.append(category_percent)

        # Append totals for each category to end
        for sample_category, sample_category_total in self._sample_categories_and_totals:
            data.append(sample_category_total)

        return data

    def summary_names(self) -> List[str]:
        return self._summary_names


class FeaturesFromIndexComparator(FeaturesComparator, abc.ABC):

    def __init__(self, connection: DataIndexConnection, include_unknown_samples: bool,
                 include_unknown_no_present_samples: bool):
        super().__init__(connection=connection, include_unknown_samples=include_unknown_samples,
                         include_unknown_no_present_samples=include_unknown_no_present_samples)

    def _get_samples_or_empty(self, feature: QueryFeature,
                              feature_id_set_dict: Dict[str, SampleSet]) -> SampleSet:
        if feature.id in feature_id_set_dict:
            return feature_id_set_dict[feature.id]
        else:
            return SampleSet.create_empty()

    def _create_summary_df(self, present_features: Dict[str, FeatureSamples],
                           present_samples: SampleSet,
                           feature_samples_summarizer: FeatureSamplesSummarizer) -> pd.DataFrame:
        data = []
        total = self._get_total(present_samples)
        for feature_id in present_features:
            feature = present_features[feature_id]

            samples_in_feature = present_samples.intersection(feature.sample_ids)

            if self._include_unknown_samples and not feature.is_unknown:
                query_feature = QueryFeatureFactory.instance().create_feature(feature.query_id)
                unknown_samples_dict = self._connection.sample_service.find_unknown_sample_sets_by_features(
                    [query_feature])
                samples_unknown_in_feature = present_samples.intersection(
                    self._get_samples_or_empty(query_feature, unknown_samples_dict))
                sample_unknown_count = len(samples_unknown_in_feature)

                samples_in_feature = FeaturesFromIndexComparator.remove_unknowns_from_present_samples(
                    sample_set=samples_in_feature,
                    unknown_set=samples_unknown_in_feature,
                    sample_service=self._connection.sample_service)
            else:
                samples_unknown_in_feature = None
                sample_unknown_count = None

            sample_count = len(samples_in_feature)

            if sample_count > 0 or (self._include_unknown_no_present_samples and
                                    sample_unknown_count is not None and sample_unknown_count > 0):
                data.append(self._create_feature_sample_count_row(feature_id,
                                                                  feature=feature,
                                                                  feature_samples=samples_in_feature,
                                                                  feature_samples_unknown=samples_unknown_in_feature,
                                                                  total=total,
                                                                  feature_samples_summarizer=feature_samples_summarizer)
                            )
        summary_df = pd.DataFrame(data,
                                  columns=[self.index_name] \
                                          + self.feature_id_columns + feature_samples_summarizer.summary_names())
        return summary_df.set_index(self.index_name)

    def summary(self, sample_set: SampleSet) -> pd.DataFrame:
        samples_summarizer = FeatureSamplesSingleCategorySummarizer()
        return self._do_summary(sample_set=sample_set, feature_samples_summarizer=samples_summarizer)

    def features_comparison(self, selected_samples: SampleSet,
                            sample_categories: List[SampleSet],
                            category_prefixes: List[str] = None,
                            category_samples_threshold: int = None,
                            unit: str = 'percent') -> pd.DataFrame:
        sample_categories_in_selected = [c.intersection(selected_samples) for c in sample_categories]

        if category_samples_threshold is not None:
            filtered_categories = []
            filtered_category_prefixes = []
            for sample_category, prefix in zip(sample_categories_in_selected, category_prefixes):
                if len(sample_category) >= category_samples_threshold:
                    filtered_categories.append(sample_category)
                    filtered_category_prefixes.append(prefix)
        else:
            filtered_categories = sample_categories_in_selected
            filtered_category_prefixes = category_prefixes

        samples_summarizer = FeatureSamplesMultipleCategorySummarizer(sample_categories=filtered_categories,
                                                                      category_prefixes=filtered_category_prefixes,
                                                                      compare_kind=unit)
        return self._do_summary(sample_set=selected_samples, feature_samples_summarizer=samples_summarizer)

    @abc.abstractmethod
    def _do_summary(self, sample_set: SampleSet, feature_samples_summarizer: FeatureSamplesSummarizer) -> pd.DataFrame:
        pass

    @abc.abstractmethod
    def _create_feature_sample_count_row(self, feature_id: str, feature: FeatureSamples,
                                         feature_samples: SampleSet,
                                         feature_samples_unknown: Optional[SampleSet],
                                         total: int,
                                         feature_samples_summarizer: FeatureSamplesSummarizer) -> List[Any]:
        pass

    @classmethod
    def remove_unknowns_from_present_samples(cls, sample_set: SampleSet, unknown_set: SampleSet,
                                      sample_service: SampleService, extra_message: str = '') -> SampleSet:
        # Handle situations where a sample is both found and unknown
        # This should not occur (unless there's bugs) but this will warn you if it does
        common_unknown_found_set = sample_set.intersection(unknown_set)
        if len(common_unknown_found_set) > 0:
            # Some of the code here is to only print the top "max" number of sample names in the warning message
            max = 10
            common_ids = list(common_unknown_found_set)[:max]
            common_names = [s.name for s in sample_service.find_samples_by_ids(common_ids)]
            if len(common_unknown_found_set) > max:
                msg = f'names=[{", ".join(common_names)}, ...]'
            else:
                msg = f'names={common_names}'
            logger.warning(f'There are {len(common_unknown_found_set)} samples ({msg}) that are both unknown and found '
                           f'{extra_message}. Will set these samples to unknown.')
            sample_set = sample_set.minus(common_unknown_found_set)

        return sample_set
