import abc
from typing import Dict, Set, Optional

from storage.variant.MaskedGenomicRegions import MaskedGenomicRegions
from storage.variant.io.FeaturesReader import FeaturesReader


class NucleotideFeaturesReader(FeaturesReader):

    def __init__(self):
        super().__init__()
        self._genomic_masked_regions: Optional[Dict[str, MaskedGenomicRegions]] = None

    def _minimal_expected_columns(self) -> Set[str]:
        return {'SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE'}

    def get_genomic_masked_regions(self) -> Dict[str, MaskedGenomicRegions]:
        if self._genomic_masked_regions is None:
            self._genomic_masked_regions = self._read_genomic_masked_regions()
        return self._genomic_masked_regions

    @abc.abstractmethod
    def _read_genomic_masked_regions(self) -> Dict[str, MaskedGenomicRegions]:
        pass
