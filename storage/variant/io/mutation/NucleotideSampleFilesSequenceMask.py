from pathlib import Path
from typing import Optional
import logging

from storage.variant.io.mutation.NucleotideSampleFiles import NucleotideSampleFiles
from storage.variant.MaskedGenomicRegions import MaskedGenomicRegions
from storage.variant.util import parse_sequence_file

logger = logging.getLogger(__name__)


class NucleotideSampleFilesSequenceMask(NucleotideSampleFiles):

    def __init__(self, sample_name: str, vcf_file: Path, vcf_file_index: Optional[Path], sample_mask_sequence: Path,
                 tmp_dir: Path):
        super().__init__(sample_name=sample_name,
                         vcf_file=vcf_file,
                         vcf_file_index=vcf_file_index,
                         mask_bed_file=None,
                         preprocessed=False,
                         tmp_dir=tmp_dir)
        self._sample_mask_sequence = sample_mask_sequence

    def _preprocess_mask(self) -> Path:
        new_file = self._tmp_dir / f'{self.sample_name}.bed.gz'
        if new_file.exists():
            raise Exception(f'File {new_file} already exists')

        name, records = parse_sequence_file(self._sample_mask_sequence)
        logger.debug(f'Getting genomic masks from {self._sample_mask_sequence}')
        masked_regions = MaskedGenomicRegions.from_sequences(sequences=records)
        masked_regions.write(new_file)
        return new_file

    @classmethod
    def create(cls, sample_name: str, vcf_file: Path, sample_mask_sequence: Path,
               tmp_dir: Path) -> NucleotideSampleFiles:
        return NucleotideSampleFilesSequenceMask(sample_name=sample_name,
                                                 vcf_file=vcf_file,
                                                 vcf_file_index=None,
                                                 sample_mask_sequence=sample_mask_sequence,
                                                 tmp_dir=tmp_dir)
