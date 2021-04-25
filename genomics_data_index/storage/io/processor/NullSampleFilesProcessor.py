from __future__ import annotations

from typing import Generator, List

from genomics_data_index.storage.io.SampleData import SampleData
from genomics_data_index.storage.io.SampleFilesProcessor import SampleFilesProcessor


class NullSampleFilesProcessor(SampleFilesProcessor):
    files_processor_instance = None

    def __init__(self):
        super().__init__()

    def process(self, sample_data: List[SampleData]) -> Generator[SampleData, None, None]:
        yield from sample_data

    @classmethod
    def instance(cls) -> NullSampleFilesProcessor:
        if cls.files_processor_instance is None:
            cls.files_processor_instance = NullSampleFilesProcessor()
        return cls.files_processor_instance
