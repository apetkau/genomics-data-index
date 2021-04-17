import abc
from typing import Dict
from pathlib import Path
import logging

from storage.variant.io.SampleFiles import SampleFiles
from storage.variant.io.SampleFilesProcessor import SampleFilesProcessor

logger = logging.getLogger(__name__)


class SerialSampleFilesProcessor(SampleFilesProcessor):

    def __init__(self, preprocess_dir: Path):
        super().__init__()
        self._preprocess_dir = preprocess_dir

    def preprocess_files(self) -> Dict[str, SampleFiles]:
        processed_files = {}
        for sample_files in self.sample_files_list():
            logger.debug(f'Pre-processing files for sample [{sample_files.sample_name}]')
            processed_files[sample_files.sample_name] = sample_files.persist(self._preprocess_dir)

        return processed_files
