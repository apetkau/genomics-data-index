import logging
import multiprocessing as mp
from pathlib import Path
from typing import Dict

from storage.variant.io.SampleFiles import SampleFiles
from storage.variant.io.SampleFilesProcessor import SampleFilesProcessor

logger = logging.getLogger(__name__)


class MultipleProcessSampleFilesProcessor(SampleFilesProcessor):

    def __init__(self, preprocess_dir: Path, processing_cores: int, max_chunk_size=1000):
        super().__init__()
        self._preprocess_dir = preprocess_dir
        self._processing_cores = processing_cores
        self._max_chunk_size = max_chunk_size

    def handle_single_file(self, sample_files: SampleFiles) -> SampleFiles:
        return sample_files.persist(self._preprocess_dir)

    def _get_chunk_size(self):
        chunk_size = max(1, int(len(self.sample_files_list()) / self._processing_cores))
        return min(self._max_chunk_size, chunk_size)

    def preprocess_files(self) -> Dict[str, SampleFiles]:
        processed_files = {}

        logger.debug(f'Starting preprocessing {len(self.sample_files_list())} samples '
                     f'with {self._processing_cores} cores')
        chunk_size = self._get_chunk_size()
        with mp.Pool(self._processing_cores) as pool:
            output_files = pool.imap_unordered(self.handle_single_file,
                                               self.sample_files_list(),
                                               chunk_size)
            for entry in output_files:
                logger.debug(f'Finished {entry.sample_name}')
                processed_files[entry.sample_name] = entry

        logger.debug(f'Finished preprocessing {len(self.sample_files_list())} samples')

        return processed_files
