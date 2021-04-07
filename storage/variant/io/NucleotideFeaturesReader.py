import abc
from typing import Dict, Set

from storage.variant.io.FeaturesReader import FeaturesReader
from storage.variant.MaskedGenomicRegions import MaskedGenomicRegions


class NucleotideFeaturesReader(FeaturesReader):

    def __init(self):
        super().__init__()

    def _minimal_expected_columns(self) -> Set[str]:
        return {'SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE'}

    def get_genomic_masked_regions(self) -> Dict[str, MaskedGenomicRegions]:
        return self._read_genomic_masked_regions()

    @abc.abstractmethod
    def _read_genomic_masked_regions(self) -> Dict[str, MaskedGenomicRegions]:
        pass
