import abc
from typing import Generator, Set, Dict

from storage.variant.io.SampleData import SampleData


class SampleDataPackage(abc.ABC):

    def __init__(self):
        pass

    @abc.abstractmethod
    def sample_names(self) -> Set[str]:
        pass

    def process_all_data(self) -> Dict[str, SampleData]:
        processed_data = {}
        for sample_data in self.iter_sample_data():
            processed_data[sample_data.sample_name] = sample_data

        return processed_data

    @abc.abstractmethod
    def iter_sample_data(self) -> Generator[SampleData, None, None]:
        pass
