import logging
from typing import List, Generator

import pandas as pd

from genomics_data_index.storage.io.mutation.NucleotideSampleData import NucleotideSampleData
from genomics_data_index.storage.io.mutation.variants_processor.VcfVariantsTableProcessor import \
    VcfVariantsTableProcessor
from genomics_data_index.storage.util import TRACE_LEVEL

logger = logging.getLogger(__name__)


class SerialVcfVariantsTableProcessor(VcfVariantsTableProcessor):

    def __init__(self):
        super().__init__()

    def process(self, sample_data: List[NucleotideSampleData], include_masked_regions: bool) -> Generator[
        pd.DataFrame, None, None]:
        for sample_files in sample_data:
            logger.log(TRACE_LEVEL, f'Constructing features table for sample [{sample_files.sample_name}]')
            yield sample_files.read_sample_data_features(include_masked_regions=include_masked_regions)
