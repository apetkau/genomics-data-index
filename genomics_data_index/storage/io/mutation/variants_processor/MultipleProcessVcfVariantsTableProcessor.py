import logging
import multiprocessing as mp
from typing import List, Generator

import pandas as pd

from genomics_data_index.storage.io.mutation.NucleotideSampleData import NucleotideSampleData
from genomics_data_index.storage.io.mutation.variants_processor.VcfVariantsTableProcessor import \
    VcfVariantsTableProcessor, VcfVariantsTableProcessorFactory

logger = logging.getLogger(__name__)


class MultipleProcessVcfVariantsTableProcessor(VcfVariantsTableProcessor):

    def __init__(self, include_masked_regions: bool, ncores: int, max_chunk_size: int):
        super().__init__(include_masked_regions)
        self._ncores = ncores
        self._max_chunk_size = max_chunk_size

    def _get_chunk_size(self, number_samples: int):
        chunk_size = max(1, int(number_samples / self._ncores))
        return min(self._max_chunk_size, chunk_size)

    def _read_sample_data(self, sample_data: NucleotideSampleData) -> pd.DataFrame:
        return sample_data.read_sample_data_features(include_masked_regions=self.include_masked_regions)

    def process(self, sample_data: List[NucleotideSampleData]) -> Generator[pd.DataFrame, None, None]:
        number_samples = len(sample_data)
        chunk_size = self._get_chunk_size(number_samples)
        logger.debug(f'Starting to read features for {number_samples} samples '
                     f'with {self._ncores} cores and chunk size {chunk_size}')

        with mp.Pool(self._ncores) as pool:
            frame_vcf_masks = pool.imap_unordered(self._read_sample_data,
                                                  sample_data,
                                                  chunk_size)
            for frame_vcf_mask in frame_vcf_masks:
                yield frame_vcf_mask

        logger.debug(f'Finished reading features table from {number_samples} samples')


class MultipleProcessVcfVariantsTableProcessorFactory(VcfVariantsTableProcessorFactory):

    def __init__(self, ncores: int = 1, max_chunk_size: int = 500):
        super().__init__()
        self._ncores = ncores
        self._max_chunk_size = max_chunk_size

    def create_variants_processor(self, include_masked_regions: bool) -> VcfVariantsTableProcessor:
        return MultipleProcessVcfVariantsTableProcessor(include_masked_regions,
                                                        ncores=self._ncores,
                                                        max_chunk_size=self._max_chunk_size)
