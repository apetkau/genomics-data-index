import logging
import os
from pathlib import Path
from typing import List, Dict, Optional

import pandas as pd
import vcf

from genomics_data_index.storage.MaskedGenomicRegions import MaskedGenomicRegions
from genomics_data_index.storage.io.SampleData import SampleData
from genomics_data_index.storage.io.mutation.NucleotideFeaturesReader import NucleotideFeaturesReader
from genomics_data_index.storage.io.mutation.NucleotideSampleData import NucleotideSampleData

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

    def _parse_annotations(self, file: Path, reader: vcf.Reader, vcf_df: pd.DataFrame) -> pd.DataFrame:
        infos = reader.infos
        expected_start_str = "Functional annotations: '"
        if 'ANN' in infos:
            ann_info = infos['ANN'].desc
            if not ann_info.startswith(expected_start_str):
                logger.warning(
                    f"Found 'ANN' in VCF file [{file}] but INFO description does not start with [{expected_start_str}]."
                    " Will ignore ANN annotations.")
            else:
                ann_info = ann_info[len(expected_start_str):]

                if ann_info.endswith("'"):
                    ann_info = ann_info[0:len(ann_info) - 1]

                ann_headers = [x.strip() for x in ann_info.split('|')]

                vcf_df_with_keys = vcf_df.copy()
                vcf_df_with_keys['VARIANT_KEY'] = vcf_df_with_keys.apply(lambda x: f"{x['CHROM']}:{x['POS']}:{x['REF']}:{x['ALT']}",
                                                     axis='columns')
                vcf_df_with_keys = vcf_df_with_keys.reset_index().rename({'index': 'original_index'}, axis='columns')
                vcf_df_with_keys = vcf_df_with_keys.set_index('VARIANT_KEY')

                # Split up multiple entries in ANN
                ann_groups = vcf_df_with_keys['INFO'].map(lambda x: x['ANN']).explode()
                print("ANN Groups\n")
                print(ann_groups)
                ann_split_fields = ann_groups.map(lambda x: x.split('|'))
                ann_split_fields = ann_split_fields.map(lambda x: dict(zip(ann_headers, x))).to_frame()
                ann_split_fields = ann_split_fields.rename({'INFO': 'ANN'}, axis='columns').reset_index()

                def insert_key_to_dictionary(x: pd.Series):
                    """
                    I use this to insert the 'VARIANT_KEY' into the 'ANN' dictionary
                    so that when I run pd.json_normalize I can keep everything linked by the same key
                    (since json_normalize removes my index in the Series).
                    :param x: A series of columns.
                    :return: A series of columns with the 'VARIANT_KEY' inserted into the 'ANN' dictionary.
                    """
                    x['ANN']['VARIANT_KEY'] = x['VARIANT_KEY']
                    return x

                ann_series = ann_split_fields.apply(insert_key_to_dictionary, axis='columns')['ANN']
                ann_df = pd.json_normalize(ann_series).add_prefix('ANN.').set_index('ANN.VARIANT_KEY')

                # Now join ann_df back to vcf_df_with_keys to recover the original index
                vcf_df_with_keys = vcf_df_with_keys.merge(ann_df, left_index=True, right_index=True).set_index('original_index')
                print("vcf_df_with_keys\n")
                print(vcf_df_with_keys)

                return vcf_df_with_keys[
                    ['ANN.Annotation', 'ANN.Gene_Name', 'ANN.Gene_ID', 'ANN.Feature_Type', 'ANN.HGVS.c', 'ANN.HGVS.p']]

        return pd.DataFrame()

    def get_or_create_feature_file(self, sample_name: str):
        vcf_file, index_file = self._sample_files_map[sample_name].get_vcf_file()
        return vcf_file

    def read_vcf(self, file: Path, sample_name: str) -> pd.DataFrame:
        reader = vcf.Reader(filename=str(file))
        df = pd.DataFrame([vars(r) for r in reader])
        out = self._fix_df_columns(df)
        print("Main DF\n")
        print(out)

        out['ALT'] = out['ALT'].map(self._fix_alt)
        out['REF'] = out['REF'].map(self._fix_ref)
        out['TYPE'] = self._get_type(out)

        ann_df = self._parse_annotations(file=file, reader=reader, vcf_df=out)
        if len(ann_df) > 0:
            out = out.merge(ann_df, how='inner', left_index=True, right_index=True)
        else:
            logger.debug(f'No snpeff annotations found (no INFO.ANN entry) in VCF [{file}].')

        out = self._drop_extra_columns(out)

        out['FILE'] = os.path.basename(file)
        cols = out.columns.tolist()
        out['SAMPLE'] = sample_name
        out = out.reindex(columns=['SAMPLE'] + cols)
        return out.loc[:, ['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE',
                           'ANN.Annotation', 'ANN.Gene_Name', 'ANN.Gene_ID', 'ANN.Feature_Type',
                           'ANN.HGVS.c', 'ANN.HGVS.p', 'FILE']]

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
    def create(cls, sample_files_map: Dict[str, NucleotideSampleData]):
        return cls(sample_files_map=sample_files_map)
