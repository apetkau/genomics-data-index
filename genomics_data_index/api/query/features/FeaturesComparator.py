import abc
from typing import List

import pandas as pd

from genomics_data_index.configuration.connector.DataIndexConnection import DataIndexConnection
from genomics_data_index.storage.SampleSet import SampleSet


class FeaturesComparator(abc.ABC):
    FEATURES_SELECTIONS = ['all', 'unique']

    def __init__(self, connection: DataIndexConnection):
        self._connection = connection

    @property
    @abc.abstractmethod
    def feature_id_columns(self) -> List[str]:
        pass

    @property
    @abc.abstractmethod
    def summary_columns(self) -> List[str]:
        pass

    @property
    @abc.abstractmethod
    def index_name(self) -> str:
        pass

    def _join_additional_columns(self, features_df: pd.DataFrame) -> pd.DataFrame:
        return features_df

    @abc.abstractmethod
    def features_comparison(self, selected_samples: SampleSet,
                            sample_categories: List[SampleSet],
                            category_prefixes: List[str] = None,
                            unit: str = 'percent') -> pd.DataFrame:
        """
        Creates a dataframe which compares different categories of samples with each other with respect to features.

        For example, if kind=='mutations', compare_kind == 'percent' and there are two sample_categories then
        this will return dataframe like:

        | Mutation     | Total | Category1_percent | Category2_percent | Category1_total | Category2_total |
        | ref:100:A:T  | 10    | 50%               | 100%              | 8               | 2               |
        | ref:200:CT:C | 10    | 100%              | 0%                | 8               | 2               |
        | ...          | ...   | ...               | ...               | 8               | 2               |

        Here, "Category1_percent" is the percent of samples in Category1 that have this mutation/feature
        (50% or 4 out of 8 samples in Category1). "Category2_percent" is the percent of samples in Category2 with the
        feature (100% or 2 out of 2 samples in Category2).

        "Category1_total" and "Category2_total" are the total samples in each category. "Total" is the total
        samples in the overall query that form the universe from which we are defining "Category1" and "Category2".

        Note: since categories are defined based on sample sets, there is no enforcement that categories are
        mutually exclusive (that is, "Category1_total" + "Category2_total" will not always equal "Total"). This
        is done on purpose in case the categories you wish to compare are not mutually exclusive.

        :param selected_samples: The set of selected samples of which sample_categories will form subsets of.
        :param sample_categories: The different categories to compare.
        :param category_prefixes: The prefixes to use for the different categories (defaults to Category1, Category2, ...).
        :param unit: The type of data to compare in each category (either 'percent', 'proportion', or 'count').
        :return: A dataframe comparing each category with respect to the differences in features.
        """
        pass

    @abc.abstractmethod
    def summary(self, sample_set: SampleSet) -> pd.DataFrame:
        """
        Given a samples set, summarizes the features for all samples in this set.
        :param sample_set: The set of samples to summarize features in.
        :return: A dataframe summarizing features in this set of samples.
        """
        pass

    def unique_summary(self, sample_set: SampleSet, other_set: SampleSet) -> pd.DataFrame:
        """
        Given a samples set, summarizes the features for all samples in this set that are not in another set
        (i.e., are unique to the given sample set compared to the other set).
        :param sample_set: The set of samples to summarize features in.
        :param other_set: The set of samples features should not appear in.
        :return: A dataframe summarizing unique features in this set of samples.
        """
        features_df = self.summary(sample_set)
        features_complement_df = self.summary(other_set)
        features_merged_df = features_df.merge(features_complement_df, left_index=True, right_index=True,
                                               how='left', indicator=True, suffixes=('_x', '_y'))
        rename_dict = {col + '_x': col for col in self.summary_columns}
        features_unique_df = features_merged_df[features_merged_df['_merge'] == 'left_only'].rename(rename_dict,
                                                                                                    axis='columns')
        features_unique_df = features_unique_df[self.summary_columns]
        return self._join_additional_columns(features_unique_df)
