import logging
from typing import List, Dict, Optional

import pandas as pd

from genomics_data_index.storage.MaskedGenomicRegions import MaskedGenomicRegions
from genomics_data_index.storage.io.SampleData import SampleData
from genomics_data_index.storage.io.mutation.NucleotideFeaturesReader import NucleotideFeaturesReader
from genomics_data_index.storage.io.mutation.NucleotideSampleData import NucleotideSampleData
from genomics_data_index.storage.io.mutation.VcfSnpEffAnnotationParser import VcfSnpEffAnnotationParser
from genomics_data_index.storage.io.mutation.variants_processor.VcfVariantsTableProcessor import \
    VcfVariantsTableProcessor

logger = logging.getLogger(__name__)


class VcfVariantsReader(NucleotideFeaturesReader):
    VCF_FRAME_COLUMNS = ['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'FILE', 'VARIANT_ID']

    def __init__(self, sample_files_map: Dict[str, NucleotideSampleData], variants_processor: VcfVariantsTableProcessor,
                 include_masked_regions: bool = True):
        super().__init__()
        self._sample_files_map = sample_files_map
        self._snpeff_parser = VcfSnpEffAnnotationParser()
        self._include_masked_regions = include_masked_regions
        self._variants_processor = variants_processor

    def get_or_create_feature_file(self, sample_name: str):
        vcf_file, index_file = self._sample_files_map[sample_name].get_vcf_file()
        return vcf_file

    def _read_features_table(self) -> pd.DataFrame:
        frames = []
        sample_data_list = list(self._sample_files_map.values())
        logger.debug(f'Starting to read features table from {len(self._sample_files_map)} VCF files')
        for frame_vcf_mask in self._variants_processor.process(sample_data_list, self._include_masked_regions):
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
    def create(cls, sample_files_map: Dict[str, NucleotideSampleData], variants_processor: VcfVariantsTableProcessor,
               include_masked_regions: bool = True):
        return cls(sample_files_map=sample_files_map, variants_processor=variants_processor,
                   include_masked_regions=include_masked_regions)
