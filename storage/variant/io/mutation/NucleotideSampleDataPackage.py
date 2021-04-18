from __future__ import annotations

from pathlib import Path
from typing import List, Generator, Dict, Set

from storage.variant.io.SampleData import SampleData
from storage.variant.io.SampleDataPackage import SampleDataPackage
from storage.variant.io.SampleFilesProcessor import SampleFilesProcessor
from storage.variant.io.mutation.NucleotideSampleDataSequenceMask import NucleotideSampleDataSequenceMask
from storage.variant.io.processor.NullSampleFilesProcessor import NullSampleFilesProcessor


class NucleotideSampleDataPackage(SampleDataPackage):

    def __init__(self, sample_data: List[SampleData],
                 sample_names: Set[str],
                 sample_files_processor: SampleFilesProcessor):
        super().__init__()
        self._sample_data = sample_data
        self._sample_files_processor = sample_files_processor
        self._sample_names = sample_names

    def sample_names(self) -> Set[str]:
        return self._sample_names

    def iter_sample_data(self) -> Generator[SampleData, None, None]:
        return self._sample_files_processor.process(self._sample_data)

    @classmethod
    def create_from_sequence_masks(cls, sample_vcf_map: Dict[str, Path],
                                   masked_genomic_files_map: Dict[str, Path] = None,
                                   sample_files_processor: SampleFilesProcessor = NullSampleFilesProcessor.instance()) -> NucleotideSampleDataPackage:
        if masked_genomic_files_map is None:
            masked_genomic_files_map = {}

        sample_names_set = set()
        sample_data_list = []
        for sample_name in sample_vcf_map:
            vcf_file = sample_vcf_map[sample_name]
            if sample_name in masked_genomic_files_map:
                mask_file = masked_genomic_files_map[sample_name]
            else:
                mask_file = None

            sample_data = NucleotideSampleDataSequenceMask.create(
                sample_name=sample_name,
                vcf_file=vcf_file,
                sample_mask_sequence=mask_file
            )

            sample_data_list.append(sample_data)
            sample_names_set.add(sample_name)

        return NucleotideSampleDataPackage(sample_data=sample_data_list,
                                           sample_names=sample_names_set,
                                           sample_files_processor=sample_files_processor
                                           )

    @classmethod
    def create_from_snippy(cls, sample_dirs: List[Path],
                           sample_files_processor: SampleFilesProcessor = NullSampleFilesProcessor.instance()) -> NucleotideSampleDataPackage:

        sample_names_set = set()
        sample_data_list = []
        for d in sample_dirs:
            sample_name = d.name
            sample_data = NucleotideSampleDataSequenceMask.create(
                sample_name=sample_name,
                vcf_file=Path(d, 'snps.vcf.gz'),
                sample_mask_sequence=Path(d, 'snps.aligned.fa')
            )
            sample_data_list.append(sample_data)
            sample_names_set.add(sample_name)

        return NucleotideSampleDataPackage(sample_data=sample_data_list,
                                           sample_names=sample_names_set,
                                           sample_files_processor=sample_files_processor
                                           )
