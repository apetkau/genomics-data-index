from typing import List, Dict
import abc
from pathlib import Path
import os
import vcf
from Bio import SeqIO
import pandas as pd
import logging

from storage.variant.CoreBitMask import CoreBitMask


logger = logging.getLogger(__name__)


class VariantsReader(abc.ABC):

    def __init(self):
        pass

    @abc.abstractmethod
    def get_variants_table(self) -> pd.DataFrame:
        pass

    @abc.abstractmethod
    def get_core_masks(self) -> Dict[str, Dict[str, CoreBitMask]]:
        pass


class SnippyVariantsReader(VariantsReader):

    def __init__(self, snippy_outputs: List[Path]):
        super().__init__()
        self._snippy_outputs = snippy_outputs
        self._sample_dirs = [sample_dir for sample_dir in os.listdir(str(self._snippy_outputs))]

    def get_variants_table(self) -> pd.DataFrame:
        vcfs = [Path(d, 'snps.vcf.gz') for d in self._sample_dirs]
        frames = [self.read_vcf(f) for f in vcfs]
        return pd.concat(frames)

    def read_vcf(self, file: Path) -> pd.DataFrame:
        reader = vcf.Reader(filename=str(file))
        df = pd.DataFrame([vars(r) for r in reader])
        out = df.merge(pd.DataFrame(df.INFO.tolist()),
                       left_index=True, right_index=True)
        out = out[['CHROM', 'POS', 'REF', 'ALT', 'DP', 'QUAL', 'RO', 'AO', 'INFO']]
        out['TYPE'] = out['INFO'].map(lambda x: x['TYPE'][0])
        out = out.drop('INFO', axis='columns')
        out['ALT'] = out['ALT'].map(lambda x: str(x[0]))
        out['REF'] = out['REF'].map(lambda x: str(x[0]))
        out['AO'] = out['AO'].map(lambda x: x[0])
        cols = out.columns.tolist()
        out['FILE'] = os.path.basename(file)
        out = out.reindex(columns=['FILE'] + cols)
        return out

    def get_core_masks(self) -> Dict[str, Dict[str, CoreBitMask]]:
        aligned_fastas = [Path(d, 'snps.aligned.fa') for d in self._sample_dirs]
        core_masks = {}
        for file in aligned_fastas:
            sample_name = os.path.dirname(file)
            logger.debug(f'Loading core masks for sample=[{sample_name}]')
            core_masks[sample_name] = self.read_core_masks(Path(file))

        return core_masks

    def read_core_masks(self, file: Path) -> Dict[str, CoreBitMask]:
        sequence_masks = {}
        for record in SeqIO.parse(file, 'fasta'):
            if record.id not in sequence_masks:
                sequence_masks[record.id] = CoreBitMask(sequence=record.seq)

        return sequence_masks
