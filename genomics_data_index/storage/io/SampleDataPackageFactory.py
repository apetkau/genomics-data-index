import abc
from pathlib import Path

from genomics_data_index.storage.io.SampleDataPackage import SampleDataPackage


class SampleDataPackageFactory:

    def __init__(self):
        pass

    @abc.abstractmethod
    def create_data_package(self) -> SampleDataPackage:
        pass
