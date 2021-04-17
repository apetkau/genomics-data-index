import logging
import os
from pathlib import Path
from typing import List, Dict
import tempfile

import pandas as pd

from storage.variant.io.mutation.VcfVariantsReader import VcfVariantsReader
from storage.variant.io.mutation.NucleotideSampleFilesSequenceMask import NucleotideSampleFilesSequenceMask
from storage.variant.io.mutation.NucleotideSampleFiles import NucleotideSampleFiles

logger = logging.getLogger(__name__)


class SnippyVariantsReader(VcfVariantsReader):

    def __init__(self, sample_files_map: Dict[str, NucleotideSampleFiles]):
        super().__init__(sample_files_map=sample_files_map)

    def _get_type(self, vcf_df: pd.DataFrame) -> pd.Series:
        return vcf_df['INFO'].map(lambda x: x['TYPE'][0])

    def _fix_df_columns(self, vcf_df: pd.DataFrame) -> pd.DataFrame:
        # If no data, I still want certain column names so that rest of code still works
        if len(vcf_df) == 0:
            vcf_df = pd.DataFrame(columns=['CHROM', 'POS', 'REF', 'ALT', 'INFO'])

        out = vcf_df.merge(pd.DataFrame(vcf_df.INFO.tolist()),
                           left_index=True, right_index=True)
        return out[['CHROM', 'POS', 'REF', 'ALT', 'INFO']]

    def _drop_extra_columns(self, vcf_df: pd.DataFrame) -> pd.DataFrame:
        return vcf_df.drop('INFO', axis='columns')


    @classmethod
    def create(cls, sample_dirs: List[Path]):
        sample_files_map = {}
        tmp_dir = Path(tempfile.mkdtemp())

        for d in sample_dirs:
            sample_name = d.name
            sample_files_map[sample_name] = NucleotideSampleFilesSequenceMask.create(
                sample_name=sample_name,
                vcf_file=Path(d, 'snps.vcf.gz'),
                sample_mask_sequence=Path(d, 'snps.aligned.fa')
            ).persist(tmp_dir)

        return SnippyVariantsReader(sample_files_map=sample_files_map)
