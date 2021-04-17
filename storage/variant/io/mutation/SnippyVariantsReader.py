from __future__ import annotations

import logging
from pathlib import Path
from typing import List, Dict

from storage.variant.io.SampleFilesProcessor import SampleFilesProcessor
from storage.variant.io.mutation.NucleotideSampleFiles import NucleotideSampleFiles
from storage.variant.io.mutation.NucleotideSampleFilesSequenceMask import NucleotideSampleFilesSequenceMask
from storage.variant.io.mutation.VcfVariantsReader import VcfVariantsReader
from storage.variant.io.processor.NullSampleFilesProcessor import NullSampleFilesProcessor

logger = logging.getLogger(__name__)


class SnippyVariantsReader(VcfVariantsReader):

    def __init__(self, sample_files_map: Dict[str, NucleotideSampleFiles]):
        super().__init__(sample_files_map=sample_files_map)

    @classmethod
    def create(cls, sample_dirs: List[Path],
               sample_files_processor: SampleFilesProcessor = NullSampleFilesProcessor()) -> SnippyVariantsReader:

        for d in sample_dirs:
            sample_name = d.name
            sample_files = NucleotideSampleFilesSequenceMask.create(
                sample_name=sample_name,
                vcf_file=Path(d, 'snps.vcf.gz'),
                sample_mask_sequence=Path(d, 'snps.aligned.fa')
            )

            sample_files_processor.add(sample_files)

        sample_files_map = sample_files_processor.preprocess_files()
        return SnippyVariantsReader(sample_files_map=sample_files_map)
