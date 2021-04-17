import abc
from typing import Dict, List

from storage.variant.io.SampleFiles import SampleFiles


class SampleFilesProcessor(abc.ABC):

    def __init__(self):
        self._sample_files_list = []

    def add(self, sample_files: SampleFiles) -> None:
        self._sample_files_list.append(sample_files)

    def sample_files_list(self) -> List[SampleFiles]:
        return self._sample_files_list

    @abc.abstractmethod
    def preprocess_files(self) -> Dict[str, SampleFiles]:
        pass
