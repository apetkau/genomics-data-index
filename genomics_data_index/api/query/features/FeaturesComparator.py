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
                            category_names: List[str] = None,
                            compare_kind: str = 'percent') -> pd.DataFrame:
        """
        Creates a dataframe which compares different categories of samples with each other with respect to features.

        For example, if compare_kind == 'percent' and there are two sample_categories then
        this will return dataframe like:

        | Mutation     | Total | Category1 | Category2 |
        | ref:100:A:T  | 10    | 50%       | 70%       |
        | ref:200:CT:C | 10    | 100%      | 0%        |
        | ...          | ...   | ...       | ...       |

        :param selected_samples: The set of selected samples of which sample_categories will form subsets of.
        :param sample_categories: The different categories to compare.
        :param category_names: The names to use for the different categories (defaults to Category1, Category2, ...).
        :param compare_kind: The type of data to compare in each category (either 'percent', or 'count').
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