import abc
from typing import Union

from genomics_data_index.api.query.SamplesQuery import SamplesQuery
from genomics_data_index.storage.SampleSet import SampleSet


class ClusterScoreMethod(abc.ABC):

    def __init__(self):
        pass

    @abc.abstractmethod
    def score(self, data: Union[SamplesQuery, SampleSet], universe_samples: SamplesQuery) -> float:
        pass
