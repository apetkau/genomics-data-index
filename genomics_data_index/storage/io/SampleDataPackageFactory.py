import abc
from typing import Generator

from genomics_data_index.storage.io.SampleDataPackage import SampleDataPackage


class SampleDataPackageFactory:

    def __init__(self):
        pass

    @abc.abstractmethod
    def create_data_package(self) -> SampleDataPackage:
        pass

    @abc.abstractmethod
    def create_data_package_iter(self, batch_size: int = 100) -> Generator[SampleDataPackage, None, None]:
        pass
