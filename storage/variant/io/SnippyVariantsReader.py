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

        sample_files_map = {}
        for d in self._sample_dirs:
            sample_name = os.path.basename(d)
            vcf_file = Path(d, 'snps.vcf.gz')
            core_mask_file = Path(d, 'snps.aligned.fa')

            sample_files_map[sample_name] = {
                'vcf_file': vcf_file,
                'core_mask_file': core_mask_file,
            }

        super().__init__(sample_files_map=sample_files_map)
