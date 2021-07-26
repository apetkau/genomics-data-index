from typing import List, Generator

import pandas as pd

from genomics_data_index.storage.io.mutation.NucleotideSampleData import NucleotideSampleData
from genomics_data_index.storage.io.mutation.variants_processor.VcfVariantsTableProcessor import \
    VcfVariantsTableProcessor


class ServialVcfVariantsTableProcessor(VcfVariantsTableProcessor):

    def __init__(self):
        pass

    def process(self, sample_data: List[NucleotideSampleData]) -> Generator[pd.DataFrame, None, None]:
        pass
