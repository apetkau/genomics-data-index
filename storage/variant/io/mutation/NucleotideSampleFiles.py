from typing import Tuple, Optional
from pathlib import Path
import shutil
import logging

from storage.variant.io.SampleFiles import SampleFiles
from storage.variant.MaskedGenomicRegions import MaskedGenomicRegions
from storage.variant.io.mutation.VariationFile import VariationFile

logger = logging.getLogger(__name__)


class NucleotideSampleFiles(SampleFiles):

    def __init__(self, sample_name: str, vcf_file: Path, vcf_file_index: Path,
                 mask_bed_file: Optional[Path],
                 preprocessed: bool):
        super().__init__(sample_name=sample_name)
        self._vcf_file = vcf_file
        self._vcf_file_index = vcf_file_index
        self._preprocessed = preprocessed
        self._mask_bed_file = mask_bed_file

    def _preprocess_mask(self, output_dir: Path) -> Path:
        if self._mask_bed_file is None:
            mask_file = output_dir / f'{self.sample_name}.bed.gz'
            self._assert_file_not_exists(mask_file, 'Cannot preprocess data')
            mask = MaskedGenomicRegions.empty_mask()
            mask.write(mask_file)
            self._mask_bed_file = mask_file
        return self._mask_bed_file

    def _preprocess_vcf(self, output_dir: Path) -> Tuple[Path, Path]:
        new_file = output_dir / f'{self.sample_name}.vcf.gz'
        self._assert_file_not_exists(new_file, 'Cannot preprocess data')
        return VariationFile(self._vcf_file).write(new_file)

    def _do_preprocess_and_persist(self, output_dir: Path) -> SampleFiles:
        processed_vcf, processed_vcf_index = self._preprocess_vcf(output_dir)
        processed_mask = self._preprocess_mask(output_dir)

        return NucleotideSampleFiles(sample_name=self.sample_name,
                                     vcf_file=processed_vcf,
                                     vcf_file_index=processed_vcf_index,
                                     mask_bed_file=processed_mask,
                                     preprocessed=True)

    def _do_persist(self, output_dir: Path) -> SampleFiles:
        if self._mask_bed_file is None:
            raise Exception('mask_bed_file is None')

        new_vcf_file = output_dir / self._vcf_file.name
        new_vcf_index = output_dir / self._vcf_file_index.name
        new_mask_bed_file = output_dir / self._mask_bed_file.name

        self._assert_file_not_exists(new_vcf_file, 'Cannot persist data')
        self._assert_file_not_exists(new_vcf_index, 'Cannot persist data')
        self._assert_file_not_exists(new_mask_bed_file, 'Cannot persist data')

        logger.debug(f'Copying VCF and BED files to [{output_dir}] for sample [{self.sample_name}]')

        shutil.copy(self._vcf_file, new_vcf_file)
        shutil.copy(self._vcf_file_index, new_vcf_index)
        shutil.copy(self._mask_bed_file, new_mask_bed_file)

        return NucleotideSampleFiles(sample_name=self.sample_name,
                                     vcf_file=new_vcf_file,
                                     vcf_file_index=new_vcf_index,
                                     mask_bed_file=new_mask_bed_file,
                                     preprocessed=True)

    def is_preprocessed(self) -> bool:
        return self._preprocessed

    def _assert_file_not_exists(self, file: Path, error_prefix: str):
        if file.exists():
            raise Exception(f'{error_prefix} for sample [{self.sample_name}]: file [{file}] exists')

    def get_vcf_file(self, ignore_preprocessed=False) -> Tuple[Path, Path]:
        if ignore_preprocessed or self._preprocessed:
            return self._vcf_file, self._vcf_file_index
        else:
            raise Exception(f'VCF file for sample [{self.sample_name}] is not preprocessed: {self._vcf_file}')

    def get_mask(self) -> MaskedGenomicRegions:
        return MaskedGenomicRegions.from_file(self.get_mask_file())

    def get_mask_file(self) -> Path:
        if self._preprocessed and self._mask_bed_file is not None:
            return self._mask_bed_file
        else:
            raise Exception(f'Sample mask file is not preprocessed for sample [{self.sample_name}]')
