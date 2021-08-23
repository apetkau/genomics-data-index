import logging
from pathlib import Path
from typing import Optional

from genomics_data_index.storage.MaskedGenomicRegions import MaskedGenomicRegions
from genomics_data_index.storage.io.mutation.NucleotideSampleData import NucleotideSampleData
from genomics_data_index.storage.util import TRACE_LEVEL

logger = logging.getLogger(__name__)


class NucleotideSampleDataBedMask(NucleotideSampleData):

    def __init__(self, sample_name: str, vcf_file: Path, vcf_file_index: Optional[Path],
                 mask_bed_file: Path):
        super().__init__(sample_name=sample_name,
                         vcf_file=vcf_file,
                         vcf_file_index=vcf_file_index,
                         mask_bed_file=mask_bed_file,
                         preprocessed=False)

    def _preprocess_mask(self, output_dir: Path) -> Path:
        if self._mask_bed_file is None:
            mask_file = output_dir / f'{self.sample_name_persistence}.bed.gz'
            mask = MaskedGenomicRegions.empty_mask()
            mask.write(mask_file)
            return mask_file
        else:
            mask_file = output_dir / f'{self.sample_name_persistence}.bed.gz'
            if mask_file.exists():
                raise Exception(f'File {mask_file} already exists')

            logger.log(TRACE_LEVEL, f'Saving genomic masks for {self._mask_bed_file}')
            persisted_mask_file = MaskedGenomicRegions.from_file(file=self._mask_bed_file)
            persisted_mask_file.write(mask_file)
            return mask_file

    @classmethod
    def create(cls, sample_name: str, vcf_file: Path, mask_bed_file: Optional[Path]) -> NucleotideSampleData:
        return NucleotideSampleDataBedMask(sample_name=sample_name, vcf_file=vcf_file,
                                           vcf_file_index=None, mask_bed_file=mask_bed_file)
