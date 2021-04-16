import logging
import os
from pathlib import Path
from typing import List, Dict

import pandas as pd
import vcf

from storage.variant.MaskedGenomicRegions import MaskedGenomicRegions
from storage.variant.io.mutation.NucleotideFeaturesReader import NucleotideFeaturesReader
from storage.variant.util import parse_sequence_file

logger = logging.getLogger(__name__)


class VcfVariantsReader(NucleotideFeaturesReader):

    def __init__(self, sample_vcf_map: Dict[str, Path],
                 masked_genomic_files_map: Dict[str, Path] = None):
        super().__init__()

        if masked_genomic_files_map is None:
            masked_genomic_files_map: Dict[str, Path] = {}

        self._sample_vcf_map = sample_vcf_map
        self._genomic_mask_files_map = masked_genomic_files_map

    def _fix_df_columns(self, vcf_df: pd.DataFrame) -> pd.DataFrame:
        # If no data, I still want certain column names so that rest of code still works
        if len(vcf_df) == 0:
            vcf_df = pd.DataFrame(columns=['CHROM', 'POS', 'REF', 'ALT'])

        return vcf_df.loc[:, ['CHROM', 'POS', 'REF', 'ALT']]

    def _drop_extra_columns(self, vcf_df: pd.DataFrame) -> pd.DataFrame:
        return vcf_df

    def get_or_create_feature_file(self, sample_name: str):
        return self._sample_vcf_map[sample_name]

    def sample_feature_files(self) -> Dict[str, Path]:
        return self._sample_vcf_map

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

        if len(vcf_df) > 0:
            return vcf_df.apply(lambda x: select_type(x['REF'], x['ALT']), axis='columns')
        else:
            return pd.Series(dtype='object')

    def _read_features_table(self) -> pd.DataFrame:
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

    def get_genomic_masked_region(self, sample_name: str) -> MaskedGenomicRegions:
        if sample_name in self._genomic_mask_files_map:
            return self.read_genomic_masks_from_file(self._genomic_mask_files_map[sample_name])
        else:
            return MaskedGenomicRegions.empty_mask()

    def read_genomic_masks_from_file(self, file: Path) -> MaskedGenomicRegions:
        name, records = parse_sequence_file(file)
        logger.debug(f'Getting genomic masks from {file}')
        return MaskedGenomicRegions.from_sequences(sequences=records)

    def samples_list(self) -> List[str]:
        return list(self._sample_vcf_map.keys())
