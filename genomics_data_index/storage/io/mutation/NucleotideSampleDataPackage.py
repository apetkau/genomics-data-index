from __future__ import annotations

from pathlib import Path
from typing import List, Generator, Dict, Set

from genomics_data_index.storage.io.SampleData import SampleData
from genomics_data_index.storage.io.SampleDataPackage import SampleDataPackage
from genomics_data_index.storage.io.SampleFilesProcessor import SampleFilesProcessor
from genomics_data_index.storage.io.mutation.NucleotideSampleDataSequenceMask import NucleotideSampleDataSequenceMask
from genomics_data_index.storage.io.processor.NullSampleFilesProcessor import NullSampleFilesProcessor


class NucleotideSampleDataPackage(SampleDataPackage):

    def __init__(self, sample_data_dict: Dict[str, SampleData],
                 sample_files_processor: SampleFilesProcessor,
                 index_unknown_missing: bool = True):
        super().__init__(index_unknown_missing=index_unknown_missing)
        self._sample_data_dict = sample_data_dict
        self._sample_files_processor = sample_files_processor
        self._sample_names = set(self._sample_data_dict.keys())

    def sample_names(self) -> Set[str]:
        return self._sample_names

    def iter_sample_data(self) -> Generator[SampleData, None, None]:
        return self._sample_files_processor.process(list(self._sample_data_dict.values()))

    @classmethod
    def create_from_sequence_masks(cls, sample_vcf_map: Dict[str, Path],
                                   masked_genomic_files_map: Dict[str, Path] = None,
                                   sample_files_processor: SampleFilesProcessor = NullSampleFilesProcessor.instance(),
                                   index_unknown_missing: bool = True) -> NucleotideSampleDataPackage:
        if masked_genomic_files_map is None:
            masked_genomic_files_map = {}

        sample_data_dict = {}
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

            sample_data_dict[sample_name] = sample_data

        return NucleotideSampleDataPackage(sample_data_dict=sample_data_dict,
                                           sample_files_processor=sample_files_processor,
                                           index_unknown_missing=index_unknown_missing
                                           )

    @classmethod
    def create_from_snippy(cls, sample_dirs: List[Path],
                           sample_files_processor: SampleFilesProcessor = NullSampleFilesProcessor.instance(),
                           index_unknown_missing: bool = True) -> NucleotideSampleDataPackage:

        sample_data_dict = {}
        for d in sample_dirs:
            sample_name = d.name
            sample_data = NucleotideSampleDataSequenceMask.create(
                sample_name=sample_name,
                vcf_file=Path(d, 'snps.vcf.gz'),
                sample_mask_sequence=Path(d, 'snps.aligned.fa')
            )
            sample_data_dict[sample_name] = sample_data

        return NucleotideSampleDataPackage(sample_data_dict=sample_data_dict,
                                           sample_files_processor=sample_files_processor,
                                           index_unknown_missing=index_unknown_missing
                                           )
