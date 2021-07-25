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
from genomics_data_index.storage.io.mutation.VcfSnpEffAnnotationParser import VcfSnpEffAnnotationParser
from genomics_data_index.storage.model import NUCLEOTIDE_UNKNOWN, NUCLEOTIDE_UNKNOWN_TYPE
from genomics_data_index.storage.util import TRACE_LEVEL

logger = logging.getLogger(__name__)


class VcfVariantsReader(NucleotideFeaturesReader):
    VCF_FRAME_COLUMNS = ['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'FILE', 'VARIANT_ID']

    def __init__(self, sample_files_map: Dict[str, NucleotideSampleData],
                 include_masked_regions: bool = True):
        super().__init__()
        self._sample_files_map = sample_files_map
        self._snpeff_parser = VcfSnpEffAnnotationParser()
        self._include_masked_regions = include_masked_regions

    def _fix_df_columns(self, vcf_df: pd.DataFrame) -> pd.DataFrame:
        # If no data, I still want certain column names so that rest of code still works
        if len(vcf_df) == 0:
            vcf_df = pd.DataFrame(columns=['CHROM', 'POS', 'REF', 'ALT', 'INFO'])

        return vcf_df.loc[:, ['CHROM', 'POS', 'REF', 'ALT', 'INFO']]

    def _drop_extra_columns(self, vcf_df: pd.DataFrame) -> pd.DataFrame:
        return vcf_df.drop('INFO', axis='columns')

    def get_or_create_feature_file(self, sample_name: str):
        vcf_file, index_file = self._sample_files_map[sample_name].get_vcf_file()
        return vcf_file

    def read_vcf(self, file: Path, sample_name: str) -> pd.DataFrame:
        reader = vcf.Reader(filename=str(file))
        df = pd.DataFrame([vars(r) for r in reader])
        out = self._fix_df_columns(df)

        out['ALT'] = out['ALT'].map(self._fix_alt)
        out['REF'] = out['REF'].map(self._fix_ref)
        out['TYPE'] = self._get_type(out)

        snpeff_headers = self._snpeff_parser.parse_annotation_headers(vcf_info=reader.infos)
        ann_df = self._snpeff_parser.parse_annotation_entries(vcf_ann_headers=snpeff_headers, vcf_df=out)
        out = out.merge(ann_df, how='left', left_index=True, right_index=True)

        out = self._drop_extra_columns(out)

        out['FILE'] = os.path.basename(file)
        cols = out.columns.tolist()
        out['SAMPLE'] = sample_name
        out = out.reindex(columns=['SAMPLE'] + cols)
        return out.loc[:, self.VCF_FRAME_COLUMNS + self._snpeff_parser.ANNOTATION_COLUMNS]

    def _get_type(self, vcf_df: pd.DataFrame) -> pd.Series:
        return vcf_df['INFO'].map(lambda x: x['TYPE'][0])

    def _read_sample_table(self, sample_name: str) -> pd.DataFrame:
        vcf_file, index_file = self._sample_files_map[sample_name].get_vcf_file()
        frame = self.read_vcf(vcf_file, sample_name)
        frame = self._snpeff_parser.select_variant_annotations(frame)

        if self._include_masked_regions:
            logger.log(TRACE_LEVEL, f'Creating unknown/missing features for sample=[{sample_name}]')
            frame_mask = self.mask_to_features(self._sample_files_map[sample_name].get_mask())
            frame_mask['SAMPLE'] = sample_name
            frame_mask['FILE'] = vcf_file.name
            logger.log(TRACE_LEVEL, f'Combining VCF and unknown/missing (mask) dataframes for sample=[{sample_name}]')
            frame_vcf_mask = self.combine_vcf_mask(frame, frame_mask)
        else:
            frame_vcf_mask = frame

        return frame_vcf_mask

    def _read_features_table(self) -> pd.DataFrame:
        frames = []
        logger.debug(f'Starting to read features table from {len(self._sample_files_map)} VCF files')
        for sample in self._sample_files_map:
            frame_vcf_mask = self._read_sample_table(sample)
            frames.append(frame_vcf_mask)

        logger.debug(f'Finished reading features table from {len(self._sample_files_map)} VCF files')
        return pd.concat(frames)

    def mask_to_features(self, genomic_mask: MaskedGenomicRegions) -> pd.DataFrame:
        mask_features = []
        ref = 1
        alt = NUCLEOTIDE_UNKNOWN
        type = NUCLEOTIDE_UNKNOWN_TYPE
        for sequence_name, position in genomic_mask.positions_iter(start_position_index='1'):
            variant_id = f'{sequence_name}:{position}:{ref}:{alt}'
            mask_features.append([sequence_name, position, ref, alt, type, variant_id])

        return pd.DataFrame(mask_features, columns=['CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'VARIANT_ID'])

    def combine_vcf_mask(self, vcf_frame: pd.DataFrame, mask_frame: pd.DataFrame) -> pd.DataFrame:
        """
        Combine features together for VCF variants dataframe with mask dataframe, checking for any overlaps.
        If there is an overlap (e.g., a variant call also is in a masked out region) the masked position (unknown/missing)
        will be preferred as the true feature.
        :param vcf_frame: The dataframe containing only mutations/variant calls.
        :param mask_frame: The dataframe containing features from the genome mask (missing/unknown positions.
        :return: The combined data frame of both types of features.
        """
        combined_df = pd.concat([vcf_frame, mask_frame])

        # Define an order column for TYPE so I can select NUCLEOTIDE_UNKNOWN_TYPE ahead of any other type
        combined_df['TYPE_ORDER'] = 1
        combined_df.loc[combined_df['TYPE'] == NUCLEOTIDE_UNKNOWN_TYPE, 'TYPE_ORDER'] = 0

        # For any overlapping positions, prefer the NUCLEOTIDE_UNKNOWN_TYPE type
        # This may not handle every potential case where a variant overlaps with a masked region
        # (e.g., indel veriants which impact more than one nucleotide) but those should not show up
        # in a VCF file AND also in the mask file if everything was called properly.
        combined_df = combined_df.sort_values(
            ['CHROM', 'POS', 'TYPE_ORDER']).groupby(['CHROM', 'POS'], sort=False).nth(0).reset_index()

        return combined_df.loc[:, self.VCF_FRAME_COLUMNS + self._snpeff_parser.ANNOTATION_COLUMNS]

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
    def create(cls, sample_files_map: Dict[str, NucleotideSampleData],
               include_masked_regions: bool = True):
        return cls(sample_files_map=sample_files_map, include_masked_regions=include_masked_regions)
