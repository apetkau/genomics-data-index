from __future__ import annotations
from typing import List, Generator, Dict, Set
from pathlib import Path

from storage.variant.io.SampleDataPackage import SampleDataPackage
from storage.variant.io.SampleData import SampleData
from storage.variant.io.SampleFilesProcessor import SampleFilesProcessor
from storage.variant.io.processor.NullSampleFilesProcessor import NullSampleFilesProcessor
from storage.variant.io.mutation.NucleotideSampleDataSequenceMask import NucleotideSampleDataSequenceMask


class NucleotideSampleDataPackage(SampleDataPackage):

    def __init__(self, sample_files_processor: SampleFilesProcessor):
        super().__init__()
        self._sample_files_processor = sample_files_processor

    # def add(self, sample_files: SampleData) -> None:
    #     self._sample_files_list.append(sample_files)
    #
    def sample_names(self) -> Set[str]:
        return {}

    def iter_sample_data(self) -> Generator[SampleData, None, None]:
        return self._sample_files_processor.preprocess_files()

    @classmethod
    def create_from_sequence_masks(cls, sample_vcf_map: Dict[str, Path],
                                   masked_genomic_files_map: Dict[str, Path] = None,
                                   sample_files_processor: SampleFilesProcessor = NullSampleFilesProcessor()) -> NucleotideSampleDataPackage:
        if masked_genomic_files_map is None:
            masked_genomic_files_map = {}

        sample_files_map = {}
        for sample_name in sample_vcf_map:
            vcf_file = sample_vcf_map[sample_name]
            if sample_name in masked_genomic_files_map:
                mask_file = masked_genomic_files_map[sample_name]
            else:
                mask_file = None

            sample_files = NucleotideSampleDataSequenceMask.create(
                sample_name=sample_name,
                vcf_file=vcf_file,
                sample_mask_sequence=mask_file
            )

            sample_files_map[sample_name] = sample_files
            sample_files_processor.add(sample_files)

        return NucleotideSampleDataPackage(sample_files_processor=sample_files_processor)

    @classmethod
    def create_from_snippy(cls, sample_dirs: List[Path],
                           sample_files_processor: SampleFilesProcessor = NullSampleFilesProcessor()) -> NucleotideSampleDataPackage:

        sample_files_map = {}
        for d in sample_dirs:
            sample_name = d.name
            sample_files = NucleotideSampleDataSequenceMask.create(
                sample_name=sample_name,
                vcf_file=Path(d, 'snps.vcf.gz'),
                sample_mask_sequence=Path(d, 'snps.aligned.fa')
            )
            sample_files_map[sample_name] = sample_files
            sample_files_processor.add(sample_files)

        return NucleotideSampleDataPackage(sample_files_processor=sample_files_processor)
