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
