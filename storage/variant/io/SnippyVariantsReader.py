import logging
import os
from pathlib import Path
from typing import List, Dict

import pandas as pd
import vcf
from Bio import SeqIO

from storage.variant.CoreBitMask import CoreBitMask
from storage.variant.io.VcfVariantsReader import VcfVariantsReader

logger = logging.getLogger(__name__)


class SnippyVariantsReader(VcfVariantsReader):

    def __init__(self, sample_dirs: List[Path]):
        self._sample_dirs = sample_dirs

        vcf_map = {}
        core_mask_map = {}
        for d in self._sample_dirs:
            sample_name = os.path.basename(d)
            vcf_map[sample_name] = Path(d, 'snps.vcf.gz')
            core_mask_map[sample_name] = Path(d, 'snps.aligned.fa')

        super().__init__(sample_vcf_map=vcf_map, core_mask_files_map=core_mask_map)

    def _get_type(self, vcf_df: pd.DataFrame) -> pd.Series:
        return vcf_df['INFO'].map(lambda x: x['TYPE'][0])

    def _fix_df_columns(self, vcf_df: pd.DataFrame) -> pd.DataFrame:
        out = vcf_df.merge(pd.DataFrame(vcf_df.INFO.tolist()),
                left_index=True, right_index=True)
        return out[['CHROM', 'POS', 'REF', 'ALT', 'INFO']]

    def _drop_extra_columns(self, vcf_df: pd.DataFrame) -> pd.DataFrame:
        return vcf_df.drop('INFO', axis='columns')
