from __future__ import annotations

import abc
from typing import Generator, Set

from genomics_data_index.storage.io.FeaturesReader import FeaturesReader
from genomics_data_index.storage.io.SampleData import SampleData


class SampleDataPackage(abc.ABC):

    def __init__(self, index_unknown_missing: bool):
        self._index_unknown_missing = index_unknown_missing

    @abc.abstractmethod
    def sample_names(self) -> Set[str]:
        pass

    def index_unknown_missing(self) -> bool:
        return self._index_unknown_missing

    @abc.abstractmethod
    def process_all_data(self) -> SampleDataPackage:
        pass

    @abc.abstractmethod
    def get_features_reader(self) -> FeaturesReader:
        pass

    @abc.abstractmethod
    def iter_sample_data(self) -> Generator[SampleData, None, None]:
        pass
