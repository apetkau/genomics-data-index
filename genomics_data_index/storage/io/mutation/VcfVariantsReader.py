import logging
from typing import List, Dict, Optional

import pandas as pd

from genomics_data_index.storage.MaskedGenomicRegions import MaskedGenomicRegions
from genomics_data_index.storage.io.SampleData import SampleData
from genomics_data_index.storage.io.mutation.NucleotideFeaturesReader import NucleotideFeaturesReader
from genomics_data_index.storage.io.mutation.NucleotideSampleData import NucleotideSampleData
from genomics_data_index.storage.io.mutation.VcfSnpEffAnnotationParser import VcfSnpEffAnnotationParser
from genomics_data_index.storage.io.mutation.variants_processor.VcfVariantsTableProcessor import \
    VcfVariantsTableProcessorFactory
from genomics_data_index.storage.util.SamplesProgressLogger import SamplesProgressLogger

logger = logging.getLogger(__name__)


class VcfVariantsReader(NucleotideFeaturesReader):
    VCF_FRAME_COLUMNS = ['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'FILE', 'VARIANT_ID']

    def __init__(self, sample_files_map: Dict[str, NucleotideSampleData],
                 variants_processor_factory: VcfVariantsTableProcessorFactory,
                 progress_logger: SamplesProgressLogger = None,
                 include_masked_regions: bool = True):
        super().__init__()
        self._sample_files_map = sample_files_map
        self._snpeff_parser = VcfSnpEffAnnotationParser()
        self._include_masked_regions = include_masked_regions
        self._variants_processor_factory = variants_processor_factory
        self._progress_logger = progress_logger

    def get_or_create_feature_file(self, sample_name: str):
        vcf_file, index_file = self._sample_files_map[sample_name].get_vcf_file()
        return vcf_file

    def _read_features_table(self) -> pd.DataFrame:
        frames = []
        sample_data_list = list(self._sample_files_map.values())
        num_samples = len(sample_data_list)
        logger.debug(f'Starting to read features table from {num_samples} sample files')
        variants_processor = self._variants_processor_factory.create_variants_processor(self._include_masked_regions)
        samples_processed = 0
        progress_logger = self._progress_logger
        if progress_logger is None:
            progress_logger = SamplesProgressLogger(stage_name='Features', stage_number=-1, total_samples=num_samples)
        update_every_nth_sample = progress_logger.get_update_every_nth_sample(percent_to_update=2)

        progress_logger.update_progress(0)
        for frame_vcf_mask in variants_processor.process(sample_data_list):
            frames.append(frame_vcf_mask)
            if (samples_processed % update_every_nth_sample) == 0:
                progress_logger.update_progress(samples_processed)
            samples_processed = samples_processed + 1
        progress_logger.update_progress(samples_processed)

        logger.debug(f'Finished reading features table from {num_samples} samples')
        return pd.concat(frames)

    def get_sample_files(self, sample_name: str) -> Optional[SampleData]:
        return self._sample_files_map[sample_name]

    def get_genomic_masked_region(self, sample_name: str) -> MaskedGenomicRegions:
        return self._sample_files_map[sample_name].get_mask()

    def samples_list(self) -> List[str]:
        return list(self._sample_files_map.keys())

    @classmethod
    def create(cls, sample_files_map: Dict[str, NucleotideSampleData],
               variants_processor_factory: VcfVariantsTableProcessorFactory,
               include_masked_regions: bool = True):
        return cls(sample_files_map=sample_files_map,
                   variants_processor_factory=variants_processor_factory,
                   include_masked_regions=include_masked_regions)
