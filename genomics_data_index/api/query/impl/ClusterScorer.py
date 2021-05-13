from __future__ import annotations
from typing import Union, Any, Optional, Callable

import pandas as pd

from genomics_data_index.storage.SampleSet import SampleSet
from genomics_data_index.api.query.SamplesQuery import SamplesQuery


class ClusterScorer:

    SCORE_KINDS = ['mrca_jaccard']

    def __init__(self, universe_samples: SamplesQuery):
        self._universe_samples = universe_samples

    def _score_sample_mrca_jaccard(self, data: Union[SamplesQuery, SampleSet]) -> float:
        if isinstance(data, SamplesQuery):
            data_set = data.sample_set
        elif isinstance(data, SampleSet):
            data_set = data
        else:
            raise Exception(f'Invalid type for data={data}. Got type {type(data)}. Expected {SamplesQuery.__name__}'
                            f' or {SampleSet.__name__}')

        mrca_samples = self._universe_samples.isin(data, kind='mrca')
        return data_set.jaccard_index(mrca_samples.sample_set)

    def score_samples(self, samples: Union[SamplesQuery, SampleSet], kind: str = 'mrca_jaccard') -> float:
        """
        Gives a score for how well the passed set of samples is clustered together in a tree.

        :param samples: The samples to score.
        :param kind: The kind of scoring method to use.
        :return: A score for how well the passed set of samples is clustered together in a tree.
        """
        if kind == 'mrca_jaccard':
            return self._score_sample_mrca_jaccard(samples)
        else:
            raise Exception(f'kind=[{kind}] is invalid. Must be one of {self.SCORE_KINDS}')

    def score_groupby(self, groupby_column: str, kind: str = 'mrca_jaccard',
                      min_samples_count: Optional[int] = None, max_samples_count: Optional[int] = None,
                      groupby_func: Callable[[Any], str] = None,
                      na_value: Any = None) -> pd.DataFrame:
        """
        Gives a score for how well sets of samples defined by the groupby_column cluster together in a tree.

        :param groupby_column: The column in the samples query dataframe to group by. When using this option you likely
                               will want to join your query with an external data frame defining extra metadata associated
                               with the samples.
        :param kind: The kind of scoring method to use.
        :param min_samples_count: The minimum samples in a group to include for scoring (default no minimum).
        :param max_samples_count: The maximum samples in a group to include for scoring (default no maximum).
        :param groupby_func: An optional function to apply to the groupby_column to be used to define the groups.
                             This gets passed a particular value from the groupby_column and returns a string representing
                             the group this value belongs to. For example if one value in 'groupby_column' is A.B.C then
                             groupby_func('A.B.C') => 'A' would say this particular row belongs to group 'A'.
        :param na_value: The value used to replace pd.NA with. If this is None (default) then NA values are excluded.
        :return: A dataframe containing scores and counts of samples (indexed by the groups in the groupby_column)
         for how well the passed set of samples is clustered together in a tree.
        """
        if kind not in self.SCORE_KINDS:
            raise Exception(f'kind=[{kind}] is invalid. Must be one of {self.SCORE_KINDS}')
        else:
            universe_df = self._universe_samples.toframe()

            if groupby_column not in universe_df:
                raise Exception(
                    f'column={groupby_column} not found in data frame from query={self._universe_samples}. '
                    f'Perhaps you meant to join this query with a dataframe before scoring clusters?')

            universe_sub_columns = universe_df[['Sample ID', groupby_column]]
            if na_value is not None:
                universe_sub_columns = universe_sub_columns.fillna(na_value)

            if groupby_func is None:
                groups_sample_sets_df = universe_sub_columns.groupby(groupby_column).agg(SampleSet)
            else:
                groups_sample_sets_df = universe_sub_columns.set_index(groupby_column).groupby(by=groupby_func).agg(SampleSet)
                groups_sample_sets_df.index.name = groupby_column

            groups_sample_sets_df['Sample Count'] = groups_sample_sets_df.apply(
                lambda x: len(x['Sample ID']), axis='columns')

            # Subset before scoring to not waste time scoring groups we don't need
            if min_samples_count is not None:
                groups_sample_sets_df = groups_sample_sets_df[groups_sample_sets_df['Sample Count'] >= min_samples_count]
            if max_samples_count is not None:
                groups_sample_sets_df = groups_sample_sets_df[groups_sample_sets_df['Sample Count'] <= max_samples_count]

            # Do scoring
            groups_sample_sets_df['Score'] = groups_sample_sets_df.apply(
                lambda x: self.score_samples(x['Sample ID']), axis='columns')

            return groups_sample_sets_df[['Score', 'Sample Count']]

    @classmethod
    def create(cls, universe_samples: SamplesQuery) -> ClusterScorer:
        """
        Creates a new ClusterScorer used to score sets of samples in a tree based on how well they cluster together.

        :param universe_samples: The set of samples under consideration. This should have a tree attached to it
                                 (i.e., should be of type TreeSamplesQuery).
        :return: A new ClusterScorer for the passed samples.
        """
        if not universe_samples.has_tree():
            raise Exception(f'Passed universe_samples={universe_samples} does not have an associated tree. '
                            f'Perhaps you forgot to join a tree with join_tree() or build with build_tree()?')

        return ClusterScorer(universe_samples=universe_samples)
