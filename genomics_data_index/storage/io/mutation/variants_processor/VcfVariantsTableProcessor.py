import abc
from typing import Generator, List

import pandas as pd

from genomics_data_index.storage.io.mutation.NucleotideSampleData import NucleotideSampleData


class VcfVariantsTableProcessor(abc.ABC):

    def __init__(self):
        pass

    @abc.abstractmethod
    def process(self, sample_data: List[NucleotideSampleData], include_masked_regions: bool) -> Generator[
        pd.DataFrame, None, None]:
        pass
