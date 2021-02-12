from typing import List, Dict
import logging
import os
from pathlib import Path

import pandas as pd
import vcf
from Bio import SeqIO

from storage.variant.io import VariantsReader
from storage.variant.CoreBitMask import CoreBitMask

logger = logging.getLogger(__name__)


class SnippyVariantsReader(VariantsReader):

    def __init__(self, sample_dirs: List[Path]):
        super().__init__()
        self._sample_dirs = sample_dirs

    def get_variants_table(self) -> pd.DataFrame:
        frames = []
        for directory in self._sample_dirs:
            vcf = Path(directory, 'snps.vcf.gz')
            sample_name = os.path.basename(directory)
            frame = self.read_vcf(vcf, sample_name)
            frames.append(frame)

        return pd.concat(frames)

    def _fix_alt(self, element: List[str]) -> str:
        '''
        Fix up the alternative string as the pyVCF package does not return them as a string.
        :param element: The element to fix.
        :return: The fixed element.
        '''
        return str(element[0])

    def _fix_ref(self, element: List[str]) -> str:
        '''
        Fix up the reference string as the pyVCF package does not return them as a string.
        :param element: The element to fix.
        :return: The fixed element.
        '''
        return str(element)

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

    def get_core_masks(self) -> Dict[str, Dict[str, CoreBitMask]]:
        aligned_fastas = [Path(d, 'snps.aligned.fa') for d in self._sample_dirs]
        core_masks = {}
        for file in aligned_fastas:
            sample_name = os.path.basename(os.path.dirname(file))
            logger.debug(f'Loading core masks for sample=[{sample_name}]')
            core_masks[sample_name] = self.read_core_masks(Path(file))

        return core_masks

    def read_core_masks(self, file: Path) -> Dict[str, CoreBitMask]:
        sequence_masks = {}
        for record in SeqIO.parse(file, 'fasta'):
            if record.id not in sequence_masks:
                sequence_masks[record.id] = CoreBitMask(sequence=record.seq)

        return sequence_masks
