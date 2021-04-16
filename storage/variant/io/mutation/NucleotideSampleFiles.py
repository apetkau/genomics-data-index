from abc import ABC
from typing import Tuple
from pathlib import Path
import os

from storage.variant.io.SampleFiles import SampleFiles
from storage.variant.io.mutation.VariationFile import VariationFile


class NucleotideSampleFiles(SampleFiles, ABC):

    def __init__(self, sample_name: str, vcf_file: Path, vcf_file_index: Path, preprocessed: bool,
                 tmp_dir: Path):
        super().__init__(sample_name=sample_name)
        self._vcf_file = vcf_file
        self._vcf_file_index = vcf_file_index
        self._tmp_dir = tmp_dir
        self._preprocessed = preprocessed

    def is_preprocessed(self) -> bool:
        return self._preprocessed

    def _preprocess_vcf(self) -> Tuple[Path, Path]:
        new_file = self._tmp_dir / f'{self.sample_name}.vcf.gz'
        return VariationFile(self._vcf_file).write(new_file)

    def get_vcf_file(self, ignore_preprocessed = False) -> Tuple[Path, Path]:
        if ignore_preprocessed or self._preprocessed:
            return self._vcf_file, self._vcf_file_index
        else:
            raise Exception(f'VCF file for sample [{self.sample_name}] is not preprocessed: {self._vcf_file}')
