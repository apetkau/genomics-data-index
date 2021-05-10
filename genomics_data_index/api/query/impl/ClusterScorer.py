from __future__ import annotations

from genomics_data_index.api.query import SamplesQuery
from genomics_data_index.api.query.impl.TreeSamplesQuery import TreeSamplesQuery


class ClusterScorer:

    SCORE_KINDS = ['mrca_jaccard']

    def __init__(self, universe_samples: TreeSamplesQuery):
        self._universe_samples = universe_samples

    def _score_sample_mrca_jaccard(self, data: SamplesQuery) -> float:
        mrca_samples = self._universe_samples.isin(data, kind='mrca')
        return data.sample_set.jaccard_index(mrca_samples.sample_set)

    def score_samples(self, data: SamplesQuery, kind: str = 'mrca_jaccard') -> float:
        if kind == 'mrca_jaccard':
            return self._score_sample_mrca_jaccard(data)
        else:
            raise Exception(f'kind=[{kind}] is invalid. Must be one of {self.SCORE_KINDS}')

    @classmethod
    def create(cls, universe_samples: SamplesQuery) -> ClusterScorer:
        if not isinstance(universe_samples, TreeSamplesQuery):
            raise Exception(f'Incorrect type for universe_samples. Got type={type(universe_samples)}, '
                            f'expected type={TreeSamplesQuery.__class__}')

        return ClusterScorer(universe_samples=universe_samples)
