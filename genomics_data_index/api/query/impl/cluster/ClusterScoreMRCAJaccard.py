from typing import Union

from genomics_data_index.api.query.SamplesQuery import SamplesQuery
from genomics_data_index.storage.SampleSet import SampleSet
from genomics_data_index.api.query.impl.cluster.ClusterScoreMethod import ClusterScoreMethod


class ClusterScoreMRCAJaccard(ClusterScoreMethod):

    def __init__(self):
        super().__init__()

    def score(self, data: Union[SamplesQuery, SampleSet], universe_samples: SamplesQuery) -> float:
        if isinstance(data, SamplesQuery):
            data_set = data.sample_set
        elif isinstance(data, SampleSet):
            data_set = data
        else:
            raise Exception(f'Invalid type for data={data}. Got type {type(data)}. Expected {SamplesQuery.__name__}'
                            f' or {SampleSet.__name__}')

        mrca_samples = universe_samples.isin(data, kind='mrca')
        return data_set.jaccard_index(mrca_samples.sample_set)
