import abc
from typing import Generator, Set, Dict

from genomics_data_index.storage.io.SampleData import SampleData
from genomics_data_index.storage.io.FeaturesReader import FeaturesReader


class SampleDataPackage(abc.ABC):

    def __init__(self, index_unknown_missing: bool):
        self._index_unknown_missing = index_unknown_missing

    @abc.abstractmethod
    def sample_names(self) -> Set[str]:
        pass

    def index_unknown_missing(self) -> bool:
        return self._index_unknown_missing

    def process_all_data(self) -> Dict[str, SampleData]:
        processed_data = {}
        for sample_data in self.iter_sample_data():
            processed_data[sample_data.sample_name] = sample_data

        return processed_data

    @abc.abstractmethod
    def get_features_reader(self) -> FeaturesReader:
        pass

    @abc.abstractmethod
    def iter_sample_data(self) -> Generator[SampleData, None, None]:
        pass
