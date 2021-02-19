import logging
import os
from pathlib import Path
from typing import List, Dict, Union

import pandas as pd
import vcf
from Bio import SeqIO

from storage.variant.CoreBitMask import CoreBitMask
from storage.variant.io import VariantsReader

logger = logging.getLogger(__name__)


class VcfVariantsReader(VariantsReader):

    def __init__(self, sample_vcf_map: Dict[str, Path],
                 core_mask_files_map: Dict[str, Path] = None,
                 empty_core_mask: Dict[str, CoreBitMask] = None):
        super().__init__()

        if core_mask_files_map is None and empty_core_mask is None:
            raise Exception(f'Both core_mask_files_map and empty_core_mask are None. Must set one of them.')
        elif core_mask_files_map is None:
            core_mask_files_map: Dict[str, Path] = {}

        sample_names_vcfs = set(sample_vcf_map.keys())
        sample_names_masks = set(core_mask_files_map.keys())

        if len(sample_names_masks - sample_names_vcfs) > 0:
            raise Exception(f'Missing the following sample VCFs: {sample_names_masks - sample_names_vcfs}')
        elif len(sample_names_vcfs - sample_names_masks) > 0 and empty_core_mask is None:
            raise Exception(f'Cannot have empty_core_mask unset when missing sample_names_masks. '
                            f'Missing samples: {sample_names_vcfs - sample_names_masks}')

        self._sample_vcf_map = sample_vcf_map
        self._core_mask_files_map = core_mask_files_map
        self._empty_core_mask = empty_core_mask

    def _subset_vcf_dataframe(self, vcf_df: pd.DataFrame) -> pd.DataFrame:
        return vcf_df[['CHROM', 'POS', 'REF', 'ALT']]

    def read_vcf(self, file: Path, sample_name: str) -> pd.DataFrame:
        reader = vcf.Reader(filename=str(file))
        df = pd.DataFrame([vars(r) for r in reader])
        out = df.merge(pd.DataFrame(df.INFO.tolist()),
                       left_index=True, right_index=True)

        out = self._subset_vcf_dataframe(out)

        out['ALT'] = out['ALT'].map(self._fix_alt)
        out['REF'] = out['REF'].map(self._fix_ref)
        out['TYPE'] = self._get_type(out)

        if 'INFO' in out:
            out = out.drop('INFO', axis='columns')

        out['FILE'] = os.path.basename(file)
        cols = out.columns.tolist()
        out['SAMPLE'] = sample_name
        out = out.reindex(columns=['SAMPLE'] + cols)
        return out

    def _get_type(self, vcf_df: pd.DataFrame) -> pd.Series:
        def select_type(ref: str, alt: str):
            ref = ref.upper()
            alt = alt.upper()

            if len(ref) == len(alt) and len(ref) == 1:
                return 'snp'
            elif len(ref) > len(alt) and ref.startswith(alt) or ref.endswith(alt):
                return 'del'
            elif len(alt) > len(ref) and alt.startswith(ref) or alt.endswith(ref):
                return 'ins'
            elif len(ref) > 0 and len(alt) > 0:
                return 'complex'
            else:
                raise Exception(f'Should not hit this case when defining type. ref=[{ref}], alt=[{alt}]')

        return vcf_df.apply(lambda x: select_type(x['REF'], x['ALT']), axis='columns')

    def _read_variants_table(self) -> pd.DataFrame:
        frames = []
        for sample in self._sample_vcf_map:
            vcf_file = self._sample_vcf_map[sample]
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

        for sample in self._sample_vcf_map:
            if sample in self._core_mask_files_map:
                core_masks[sample] = self.read_core_masks_from_file(self._core_mask_files_map[sample])
            else:
                core_masks[sample] = self._empty_core_mask

        return core_masks

    def read_core_masks_from_file(self, file: Path) -> Dict[str, CoreBitMask]:
        sequence_masks = {}
        for record in SeqIO.parse(file, 'fasta'):
            if record.id not in sequence_masks:
                sequence_masks[record.id] = CoreBitMask.from_sequence(sequence=record.seq)

        return sequence_masks

    def samples_list(self) -> List[str]:
        return list(self._core_mask_files_map.keys())
