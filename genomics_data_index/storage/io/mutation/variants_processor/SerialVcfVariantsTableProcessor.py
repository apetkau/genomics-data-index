from __future__ import annotations

import logging
from typing import List, Generator

import pandas as pd

from genomics_data_index.storage.io.mutation.NucleotideSampleData import NucleotideSampleData
from genomics_data_index.storage.io.mutation.variants_processor.VcfVariantsTableProcessor import \
    VcfVariantsTableProcessor, VcfVariantsTableProcessorFactory
from genomics_data_index.storage.util import TRACE_LEVEL

logger = logging.getLogger(__name__)


class SerialVcfVariantsTableProcessor(VcfVariantsTableProcessor):

    def __init__(self, include_masked_regions: bool):
        super().__init__(include_masked_regions=include_masked_regions)

    def process(self, sample_data: List[NucleotideSampleData]) -> Generator[pd.DataFrame, None, None]:
        for sample_files in sample_data:
            logger.log(TRACE_LEVEL, f'Constructing features table for sample [{sample_files.sample_name}]')
            yield sample_files.read_sample_data_features(include_masked_regions=self.include_masked_regions)


class SerialVcfVariantsTableProcessorFactory(VcfVariantsTableProcessorFactory):
    processor_factory_instance = None

    def __init__(self):
        super().__init__()

    def create_variants_processor(self, include_masked_regions: bool) -> VcfVariantsTableProcessor:
        return SerialVcfVariantsTableProcessor(include_masked_regions)

    @classmethod
    def instance(cls) -> SerialVcfVariantsTableProcessorFactory:
        if cls.processor_factory_instance is None:
            cls.processor_factory_instance = SerialVcfVariantsTableProcessorFactory()
        return cls.processor_factory_instance
