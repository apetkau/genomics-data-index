from typing import Dict
from pathlib import Path
import logging
import multiprocessing as mp

from storage.variant.io.SampleFiles import SampleFiles
from storage.variant.io.SampleFilesProcessor import SampleFilesProcessor

logger = logging.getLogger(__name__)


class MultipleProcessSampleFilesProcessor(SampleFilesProcessor):

    def __init__(self, preprocess_dir: Path, processing_cores: int):
        super().__init__()
        self._preprocess_dir = preprocess_dir
        self._processing_cores = processing_cores

    def handle_single_file(self, sample_files: SampleFiles) -> SampleFiles:
        return sample_files.persist(self._preprocess_dir)

    def preprocess_files(self) -> Dict[str, SampleFiles]:
        processed_files = {}

        with mp.Pool(self._processing_cores) as pool:
            output_files = pool.map(self.handle_single_file, self.sample_files_list())
            for entry in output_files:
                processed_files[entry.sample_name] = entry

        return processed_files
