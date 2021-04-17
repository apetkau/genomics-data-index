import abc
from typing import Generator, Set

from storage.variant.io.SampleData import SampleData


class SampleDataPackage(abc.ABC):

    def __init__(self):
        pass

    @abc.abstractmethod
    def sample_names(self) -> Set[str]:
        pass

    @abc.abstractmethod
    def iter_sample_data(self) -> Generator[SampleData, None, None]:
        pass
