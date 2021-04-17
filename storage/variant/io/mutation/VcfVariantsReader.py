import logging
import os
from pathlib import Path
from typing import List, Dict, Optional, Generator

import pandas as pd
import vcf

from storage.variant.MaskedGenomicRegions import MaskedGenomicRegions
from storage.variant.io.SampleData import SampleData
from storage.variant.io.SampleFilesProcessor import SampleFilesProcessor
from storage.variant.io.mutation.NucleotideFeaturesReader import NucleotideFeaturesReader
from storage.variant.io.mutation.NucleotideSampleData import NucleotideSampleData
from storage.variant.io.mutation.NucleotideSampleDataSequenceMask import NucleotideSampleDataSequenceMask
from storage.variant.io.processor.NullSampleFilesProcessor import NullSampleFilesProcessor
from storage.variant.io.SampleFilesProcessor import SampleFilesProcessor

logger = logging.getLogger(__name__)


class VcfVariantsReader(NucleotideFeaturesReader):

    def __init__(self, sample_files_map: Dict[str, NucleotideSampleData]):
        super().__init__()
        self._sample_files_map = sample_files_map

    def _fix_df_columns(self, vcf_df: pd.DataFrame) -> pd.DataFrame:
        # If no data, I still want certain column names so that rest of code still works
        if len(vcf_df) == 0:
            vcf_df = pd.DataFrame(columns=['CHROM', 'POS', 'REF', 'ALT', 'INFO'])

        return vcf_df.loc[:, ['CHROM', 'POS', 'REF', 'ALT', 'INFO']]

    def _drop_extra_columns(self, vcf_df: pd.DataFrame) -> pd.DataFrame:
        return vcf_df.drop('INFO', axis='columns')

    def get_or_create_feature_file(self, sample_name: str):
        vcf_file, index_file = self._sample_files_map[sample_name].get_vcf_file()
        return vcf_file

    def read_vcf(self, file: Path, sample_name: str) -> pd.DataFrame:
        reader = vcf.Reader(filename=str(file))
        df = pd.DataFrame([vars(r) for r in reader])
        out = self._fix_df_columns(df)

        out['ALT'] = out['ALT'].map(self._fix_alt)
        out['REF'] = out['REF'].map(self._fix_ref)
        out['TYPE'] = self._get_type(out)

        out = self._drop_extra_columns(out)

        out['FILE'] = os.path.basename(file)
        cols = out.columns.tolist()
        out['SAMPLE'] = sample_name
        out = out.reindex(columns=['SAMPLE'] + cols)
        return out

    def _get_type(self, vcf_df: pd.DataFrame) -> pd.Series:
        return vcf_df['INFO'].map(lambda x: x['TYPE'][0])

    def _read_features_table(self) -> pd.DataFrame:
        frames = []
        for sample in self._sample_files_map:
            vcf_file, index_file = self._sample_files_map[sample].get_vcf_file()
            frame = self.read_vcf(vcf_file, sample)
            frames.append(frame)

        return pd.concat(frames)

    def _fix_alt(self, element: List[str]) -> str:
        """
        Fix up the alternative string as the pyVCF package does not return them as a string.
        :param element: The element to fix.
        :return: The fixed element.
        """
        return str(element[0])

    def _fix_ref(self, element: List[str]) -> str:
        """
        Fix up the reference string as the pyVCF package does not return them as a string.
        :param element: The element to fix.
        :return: The fixed element.
        """
        return str(element)

    def get_sample_files(self, sample_name: str) -> Optional[SampleData]:
        return self._sample_files_map[sample_name]

    def get_genomic_masked_region(self, sample_name: str) -> MaskedGenomicRegions:
        return self._sample_files_map[sample_name].get_mask()

    def samples_list(self) -> List[str]:
        return list(self._sample_files_map.keys())

    @classmethod
    def create_from_sequence_masks(cls, sample_vcf_map: Dict[str, Path],
                                   masked_genomic_files_map: Dict[str, Path] = None,
                                   sample_files_processor: SampleFilesProcessor = NullSampleFilesProcessor()):
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

        return VcfVariantsReader(sample_files_map=sample_files_map)
