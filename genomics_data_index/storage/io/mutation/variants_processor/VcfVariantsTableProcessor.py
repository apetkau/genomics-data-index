import abc
from typing import Generator, List

import pandas as pd

from genomics_data_index.storage.io.mutation.NucleotideSampleData import NucleotideSampleData


class VcfVariantsTableProcessor(abc.ABC):

    def __init__(self, include_masked_regions: bool):
        self._include_masked_regions = include_masked_regions

    @property
    def include_masked_regions(self) -> bool:
        return self._include_masked_regions

    @abc.abstractmethod
    def process(self, sample_data: List[NucleotideSampleData]) -> Generator[pd.DataFrame, None, None]:
        pass


class VcfVariantsTableProcessorFactory(abc.ABC):

    def __init__(self):
        pass

    @abc.abstractmethod
    def create_variants_processor(self, include_masked_regions: bool) -> VcfVariantsTableProcessor:
        """
        Builds a new instance of a VcfVariantsTableProcessor.
        The reason I have a parallel factory class for each VcfVariantsTableProcessor is that for multiprocessing
        I cannot set a separate variable 'include_masked_regions' before I use Python's multiprocessing pool
        to divide up tasks. Or, I can set such a variable, but it has to be a class variable. But the decision
        on whether or not to include masked regions needs to be made much later than the decision on whether or
        not I want to use parallel/serial processing. So I create a Factory which defines the decision on whether to
        use parallel/serial processing and then the moment the decision on whether or not to use 'include_masked_regions'
        is made I can use a factory to create a new instance of VcfVariantsTableProcessor with the appropriate class
        variable set.
        :param include_masked_regions: Whether or not masked/missing/unknown regions should be included.
        :return: A new instance of a VcfVariantsTableProcessor.
        """
        pass
