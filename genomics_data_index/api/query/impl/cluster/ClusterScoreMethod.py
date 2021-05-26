import abc
from typing import Union

from genomics_data_index.api.query.SamplesQuery import SamplesQuery
from genomics_data_index.storage.SampleSet import SampleSet


class ClusterScoreMethod(abc.ABC):
    """
    Defines a method for assigning scores to clusters of samples based on their relationship in a tree.
    """

    def __init__(self):
        pass

    @abc.abstractmethod
    def score(self, data: Union[SamplesQuery, SampleSet], universe_samples: SamplesQuery) -> float:
        """
        Scores the passed set of samples (within the given universe of samples). The score is a number defining
        how well these samples cluster together in a tree. An implementation of ClusterScoreMethod can be regisered
        in :py:class:`genomics_data_index.api.query.impl.ClusterScorer` and then executed with ClusterScorer.score_samples().
        :param data: The data to score (set of samples).
        :param universe_samples: The set of universe samples.
        :return: A float defining the score of the passed set of samples in a tree.
        """
        pass
