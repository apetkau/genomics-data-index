from typing import Dict

from storage.variant.io.SampleFiles import SampleFiles
from storage.variant.io.SampleFilesProcessor import SampleFilesProcessor


class NullSampleFilesProcessor(SampleFilesProcessor):

    def __init__(self):
        super().__init__()

    def preprocess_files(self) -> Dict[str, SampleFiles]:
        processed_files = {}
        for sample_files in self.sample_files_list():
            processed_files[sample_files.sample_name] = sample_files

        return processed_files
