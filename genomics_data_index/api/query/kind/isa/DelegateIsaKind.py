import abc
from typing import Union, List

from genomics_data_index.api.query.SamplesQuery import SamplesQuery
from genomics_data_index.api.query.kind.IsaKind import IsaKind
from genomics_data_index.storage.SampleSet import SampleSet


class DelegateIsaKind(IsaKind, abc.ABC):

    def __init__(self):
        super().__init__()

    @abc.abstractmethod
    def isa(self, data: Union[str, List[str], SamplesQuery, SampleSet], query: SamplesQuery) -> SamplesQuery:
        pass
