import abc
from typing import Union, List

from genomics_data_index.api.query.SamplesQuery import SamplesQuery
from genomics_data_index.api.query.kind.isa.DelegateIsaKind import DelegateIsaKind
from genomics_data_index.storage.SampleSet import SampleSet


class SamplesTypingIsaKind(DelegateIsaKind, abc.ABC):

    def __init__(self):
        super().__init__()

    def isa(self, data: Union[str, List[str], SamplesQuery, SampleSet], query: SamplesQuery) -> SamplesQuery:
        if isinstance(data, str):
            return self.isa_type(data, query)
        else:
            raise Exception(f'Can only apply {self.__class__} to data of type string. Got data={data}')

    @abc.abstractmethod
    @property
    def name(self) -> str:
        pass

    @abc.abstractmethod
    @property
    def version(self) -> str:
        pass

    @abc.abstractmethod
    def isa_type(self, data: str, query: SamplesQuery) -> SamplesQuery:
        pass
