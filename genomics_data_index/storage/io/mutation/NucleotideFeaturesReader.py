import abc
from typing import Set

from genomics_data_index.storage.MaskedGenomicRegions import MaskedGenomicRegions
from genomics_data_index.storage.io.FeaturesReader import FeaturesReader


class NucleotideFeaturesReader(FeaturesReader):

    def __init__(self):
        super().__init__()

    def _minimal_expected_columns(self) -> Set[str]:
        return {'SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE'}

    @abc.abstractmethod
    def get_genomic_masked_region(self, sample_name: str) -> MaskedGenomicRegions:
        pass
