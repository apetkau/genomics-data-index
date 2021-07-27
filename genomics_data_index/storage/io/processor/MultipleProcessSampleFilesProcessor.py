import logging
import multiprocessing as mp
from pathlib import Path
from typing import Generator, List

from genomics_data_index.storage.io.SampleData import SampleData
from genomics_data_index.storage.io.SampleFilesProcessor import SampleFilesProcessor
from genomics_data_index.storage.util import TRACE_LEVEL

logger = logging.getLogger(__name__)


class MultipleProcessSampleFilesProcessor(SampleFilesProcessor):

    def __init__(self, preprocess_dir: Path, processing_cores: int, max_chunk_size=25):
        super().__init__()
        self._preprocess_dir = preprocess_dir
        self._processing_cores = processing_cores
        self._max_chunk_size = max_chunk_size

    def handle_single_file(self, sample_files: SampleData) -> SampleData:
        return sample_files.persist(self._preprocess_dir)

    def _get_chunk_size(self, number_samples: int):
        chunk_size = max(1, int(number_samples / self._processing_cores))
        return min(self._max_chunk_size, chunk_size)

    def process(self, sample_data: List[SampleData]) -> Generator[SampleData, None, None]:
        number_samples = len(sample_data)
        chunk_size = self._get_chunk_size(number_samples)
        logger.debug(f'Starting preprocessing {number_samples} samples '
                     f'with {self._processing_cores} cores and chunk size {chunk_size}')
        with mp.Pool(self._processing_cores) as pool:
            processed_sample_filed = pool.imap_unordered(self.handle_single_file,
                                                         sample_data,
                                                         chunk_size)
            for processed_file in processed_sample_filed:
                logger.log(TRACE_LEVEL, f'Finished {processed_file.sample_name}')
                yield processed_file

        logger.debug(f'Finished preprocessing {number_samples} samples')
