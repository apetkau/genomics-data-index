from pathlib import Path
from typing import Optional
import logging

from storage.variant.io.SampleFiles import SampleFiles
from storage.variant.io.mutation.NucleotideSampleFiles import NucleotideSampleFiles
from storage.variant.MaskedGenomicRegions import MaskedGenomicRegions
from storage.variant.util import parse_sequence_file

logger = logging.getLogger(__name__)


class NucleotideSampleFilesSequenceMask(NucleotideSampleFiles):

    def __init__(self, sample_name: str, vcf_file: Path, vcf_file_index: Optional[Path], sample_mask: Path, tmp_dir: Path,
                 preprocessed: bool = False):
        super().__init__(sample_name=sample_name, vcf_file=vcf_file, vcf_file_index=vcf_file_index,
                         preprocessed=preprocessed,
                         tmp_dir=tmp_dir)
        self._sample_mask = sample_mask

    def _preprocess_mask(self) -> Path:
        new_file = self._tmp_dir / f'{self.sample_name}.bed.gz'
        if new_file.exists():
            raise Exception(f'File {new_file} already exists')

        name, records = parse_sequence_file(self._sample_mask)
        logger.debug(f'Getting genomic masks from {self._sample_mask}')
        masked_regions = MaskedGenomicRegions.from_sequences(sequences=records)
        masked_regions.write(new_file)
        return new_file

    def _do_preprocess(self) -> SampleFiles:
        processed_vcf, processed_vcf_index = self._preprocess_vcf()
        processed_mask = self._preprocess_mask()

        return NucleotideSampleFilesSequenceMask(sample_name=self.sample_name,
                                                 vcf_file=processed_vcf,
                                                 vcf_file_index=processed_vcf_index,
                                                 sample_mask=processed_mask,
                                                 preprocessed=True,
                                                 tmp_dir=self._tmp_dir)

    @classmethod
    def create(cls, sample_name: str, vcf_file: Path, sample_mask: Path, tmp_dir: Path) -> NucleotideSampleFiles:
        return NucleotideSampleFilesSequenceMask(sample_name=sample_name,
                                                 vcf_file=vcf_file,
                                                 vcf_file_index=None,
                                                 sample_mask=sample_mask,
                                                 preprocessed=False,
                                                 tmp_dir=tmp_dir)
