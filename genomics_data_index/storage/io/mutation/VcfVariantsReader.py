import logging
import os
from pathlib import Path
from typing import List, Dict, Optional
import multiprocessing as mp

import pandas as pd
import vcf

from genomics_data_index.storage.MaskedGenomicRegions import MaskedGenomicRegions
from genomics_data_index.storage.io.SampleData import SampleData
from genomics_data_index.storage.io.mutation.NucleotideFeaturesReader import NucleotideFeaturesReader
from genomics_data_index.storage.io.mutation.NucleotideSampleData import NucleotideSampleData
from genomics_data_index.storage.io.mutation.VcfSnpEffAnnotationParser import VcfSnpEffAnnotationParser
from genomics_data_index.storage.model import NUCLEOTIDE_UNKNOWN, NUCLEOTIDE_UNKNOWN_TYPE
from genomics_data_index.storage.util import TRACE_LEVEL

logger = logging.getLogger(__name__)


class VcfVariantsReader(NucleotideFeaturesReader):
    VCF_FRAME_COLUMNS = ['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'FILE', 'VARIANT_ID']

    def __init__(self, sample_files_map: Dict[str, NucleotideSampleData], ncores: int = 1,
                 include_masked_regions: bool = True):
        super().__init__()
        self._sample_files_map = sample_files_map
        self._snpeff_parser = VcfSnpEffAnnotationParser()
        self._include_masked_regions = include_masked_regions
        self._processing_cores = ncores
        # self._variants_processor =

    def get_or_create_feature_file(self, sample_name: str):
        vcf_file, index_file = self._sample_files_map[sample_name].get_vcf_file()
        return vcf_file

    def _get_chunk_size(self, number_samples: int):
        max_chunk_size = 25
        chunk_size = max(1, int(number_samples / self._processing_cores))
        return min(max_chunk_size, chunk_size)

    def read_sample_data_features(self, sample_data: NucleotideSampleData) -> pd.DataFrame:
        return sample_data.read_sample_data_features(self._include_masked_regions)

    def _read_features_table(self) -> pd.DataFrame:
        frames = []
        number_samples = len(self._sample_files_map)
        sample_data_list = list(self._sample_files_map.values())
        chunk_size = self._get_chunk_size(number_samples)
        logger.debug(f'Starting to read features table from {len(self._sample_files_map)} VCF files')
        logger.debug(f'Starting preprocessing {number_samples} samples '
                     f'with {self._processing_cores} cores and chunk size {chunk_size}')
        with mp.Pool(self._processing_cores) as pool:
            frame_vcf_masks = pool.imap_unordered(self.read_sample_data_features,
                                                  sample_data_list,
                                                  chunk_size)
            for frame_vcf_mask in frame_vcf_masks:
                frames.append(frame_vcf_mask)

        logger.debug(f'Finished reading features table from {len(self._sample_files_map)} VCF files')
        return pd.concat(frames)

    def get_sample_files(self, sample_name: str) -> Optional[SampleData]:
        return self._sample_files_map[sample_name]

    def get_genomic_masked_region(self, sample_name: str) -> MaskedGenomicRegions:
        return self._sample_files_map[sample_name].get_mask()

    def samples_list(self) -> List[str]:
        return list(self._sample_files_map.keys())

    @classmethod
    def create(cls, sample_files_map: Dict[str, NucleotideSampleData],
               include_masked_regions: bool = True, ncores: int = 1):
        return cls(sample_files_map=sample_files_map, include_masked_regions=include_masked_regions,
                   ncores=ncores)
