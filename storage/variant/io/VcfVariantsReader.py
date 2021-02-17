import logging
import os
from pathlib import Path
from typing import List, Dict

import pandas as pd
import vcf
from Bio import SeqIO

from storage.variant.CoreBitMask import CoreBitMask
from storage.variant.io import VariantsReader

logger = logging.getLogger(__name__)


class VcfVariantsReader(VariantsReader):

    def __init__(self, sample_files_map: Dict[str, Dict[str, Path]]):
        super().__init__()
        self._sample_files_map = sample_files_map

    def read_vcf(self, file: Path, sample_name: str) -> pd.DataFrame:
        reader = vcf.Reader(filename=str(file))
        df = pd.DataFrame([vars(r) for r in reader])
        out = df.merge(pd.DataFrame(df.INFO.tolist()),
                       left_index=True, right_index=True)
        out = out[['CHROM', 'POS', 'REF', 'ALT', 'DP', 'QUAL', 'RO', 'AO', 'INFO']]
        out['TYPE'] = out['INFO'].map(lambda x: x['TYPE'][0])
        out = out.drop('INFO', axis='columns')

        out['ALT'] = out['ALT'].map(self._fix_alt)
        out['REF'] = out['REF'].map(self._fix_ref)

        out['FILE'] = os.path.basename(file)
        cols = out.columns.tolist()
        out['SAMPLE'] = sample_name
        out = out.reindex(columns=['SAMPLE'] + cols)
        return out

    def _read_variants_table(self) -> pd.DataFrame:
        frames = []
        for sample in self._sample_files_map:
            vcf_file = self._sample_files_map[sample]['vcf_file']
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

    def _read_core_masks(self) -> Dict[str, Dict[str, CoreBitMask]]:
        core_masks = {}
        for sample in self._sample_files_map:
            core_mask_file = self._sample_files_map[sample]['core_mask_file']
            core_masks[sample] = self.read_core_masks_from_file(core_mask_file)

        return core_masks

    def read_core_masks_from_file(self, file: Path) -> Dict[str, CoreBitMask]:
        sequence_masks = {}
        for record in SeqIO.parse(file, 'fasta'):
            if record.id not in sequence_masks:
                sequence_masks[record.id] = CoreBitMask.from_sequence(sequence=record.seq)

        return sequence_masks

    def samples_list(self) -> List[str]:
        return list(self._sample_files_map.keys())
