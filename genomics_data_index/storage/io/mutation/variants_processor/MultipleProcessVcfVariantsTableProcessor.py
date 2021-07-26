from typing import List, Generator

import pandas as pd

from genomics_data_index.storage.io.mutation.NucleotideSampleData import NucleotideSampleData
from genomics_data_index.storage.io.mutation.variants_processor.VcfVariantsTableProcessor import \
    VcfVariantsTableProcessor


class MultipleProcessVcfVariantsTableProcessor(VcfVariantsTableProcessor):

    def __init__(self, ncores: int = 1):
        super().__init__()
        self._ncores = ncores

    # def _get_chunk_size(self, number_samples: int):
    #     max_chunk_size = 25
    #     chunk_size = max(1, int(number_samples / self._ncores))
    #     return min(max_chunk_size, chunk_size)

    def process(self, sample_data: List[NucleotideSampleData], include_masked_regions: bool) -> Generator[
        pd.DataFrame, None, None]:
        raise NotImplementedError()

    # def _read_features_table(self) -> pd.DataFrame:
    #     frames = []
    #     number_samples = len(self._sample_files_map)
    #     sample_data_list = list(self._sample_files_map.values())
    #     chunk_size = self._get_chunk_size(number_samples)
    #     logger.debug(f'Starting to read features table from {len(self._sample_files_map)} VCF files')
    #     logger.debug(f'Starting preprocessing {number_samples} samples '
    #                  f'with {self._processing_cores} cores and chunk size {chunk_size}')
    #     with mp.Pool(self._processing_cores) as pool:
    #         frame_vcf_masks = pool.imap_unordered(self.read_sample_data_features,
    #                                               sample_data_list,
    #                                               chunk_size)
    #         for frame_vcf_mask in frame_vcf_masks:
    #             frames.append(frame_vcf_mask)
    #
    #     logger.debug(f'Finished reading features table from {len(self._sample_files_map)} VCF files')
    #     return pd.concat(frames)
