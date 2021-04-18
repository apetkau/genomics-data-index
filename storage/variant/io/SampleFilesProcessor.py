import abc
from typing import Generator, List

from storage.variant.io.SampleData import SampleData


class SampleFilesProcessor(abc.ABC):

    def __init__(self):
        pass

    @abc.abstractmethod
    def process(self, sample_data: List[SampleData]) -> Generator[SampleData, None, None]:
        pass
