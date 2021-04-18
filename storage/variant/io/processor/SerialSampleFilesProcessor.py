import logging
from pathlib import Path
from typing import Generator, List

from storage.variant.io.SampleData import SampleData
from storage.variant.io.SampleFilesProcessor import SampleFilesProcessor

logger = logging.getLogger(__name__)


class SerialSampleFilesProcessor(SampleFilesProcessor):

    def __init__(self, preprocess_dir: Path):
        super().__init__()
        self._preprocess_dir = preprocess_dir

    def process(self, sample_data: List[SampleData]) -> Generator[SampleData, None, None]:
        for sample_files in sample_data:
            logger.debug(f'Pre-processing files for sample [{sample_files.sample_name}]')
            yield sample_files.persist(self._preprocess_dir)
