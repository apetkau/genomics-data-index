import logging
from pathlib import Path
from typing import Optional

from genomics_data_index.storage.MaskedGenomicRegions import MaskedGenomicRegions
from genomics_data_index.storage.io.mutation.NucleotideSampleData import NucleotideSampleData
from genomics_data_index.storage.io.mutation.SequenceFile import SequenceFile
from genomics_data_index.storage.util import TRACE_LEVEL

logger = logging.getLogger(__name__)


class NucleotideSampleDataSequenceMask(NucleotideSampleData):

    def __init__(self, sample_name: str, vcf_file: Path, vcf_file_index: Optional[Path],
                 sample_mask_sequence: Optional[Path], subtract_vcf_from_mask: bool = False):
        super().__init__(sample_name=sample_name,
                         vcf_file=vcf_file,
                         vcf_file_index=vcf_file_index,
                         mask_bed_file=None,
                         preprocessed=False)
        self._sample_mask_sequence = sample_mask_sequence
        self._subtract_vcf_from_mask = subtract_vcf_from_mask

    def _preprocess_mask(self, output_dir: Path) -> Path:
        if self._sample_mask_sequence is None:
            mask_file = output_dir / f'{self.sample_name_persistence}.bed.gz'
            mask = MaskedGenomicRegions.empty_mask()
            mask.write(mask_file)
            return mask_file
        else:
            new_file = output_dir / f'{self.sample_name_persistence}.bed.gz'
            if new_file.exists():
                raise Exception(f'File {new_file} already exists')

            name, records = SequenceFile(self._sample_mask_sequence).parse_sequence_file()
            logger.log(TRACE_LEVEL, f'Getting genomic masks from {self._sample_mask_sequence}')
            masked_regions = MaskedGenomicRegions.from_sequences(sequences=records)

            if self._subtract_vcf_from_mask:
                logger.log(TRACE_LEVEL, f'Subtracting variants in vcf_file=[{self._vcf_file}] '
                                        f'from mask produced by sequence=[{self._sample_mask_sequence}]')
                vcf_regions = MaskedGenomicRegions.from_vcf_file(self._vcf_file)
                masked_regions = masked_regions.subtract(vcf_regions)

            masked_regions.write(new_file)
            return new_file

    @classmethod
    def create(cls, sample_name: str, vcf_file: Path,
               sample_mask_sequence: Optional[Path] = None,
               subtract_vcf_from_mask: bool = False) -> NucleotideSampleData:
        return NucleotideSampleDataSequenceMask(sample_name=sample_name,
                                                vcf_file=vcf_file,
                                                vcf_file_index=None,
                                                sample_mask_sequence=sample_mask_sequence,
                                                subtract_vcf_from_mask=subtract_vcf_from_mask)
