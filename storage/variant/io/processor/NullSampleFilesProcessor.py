from typing import Generator

from storage.variant.io.SampleFiles import SampleFiles
from storage.variant.io.SampleFilesProcessor import SampleFilesProcessor


class NullSampleFilesProcessor(SampleFilesProcessor):

    def __init__(self):
        super().__init__()

    def preprocess_files(self) -> Generator[SampleFiles, None, None]:
        yield from self.sample_files_list()
