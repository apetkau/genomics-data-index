import logging
import os
from pathlib import Path
from typing import List, Dict, Optional
import multiprocessing as mp

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

    def __init__(self, sample_files_map: Dict[str, NucleotideSampleData], ncores: int = 1,
                 include_masked_regions: bool = True):
        super().__init__()
        self._sample_files_map = sample_files_map
        self._snpeff_parser = VcfSnpEffAnnotationParser()
        self._include_masked_regions = include_masked_regions
        self._processing_cores = ncores
        # self._variants_processor =

    def get_or_create_feature_file(self, sample_name: str):
        vcf_file, index_file = self._sample_files_map[sample_name].get_vcf_file()
        return vcf_file

    def read_sample_data_features(self, sample_data: NucleotideSampleData) -> pd.DataFrame:
        vcf_file, index_file = sample_data.get_vcf_file()
        sample_name = sample_data.sample_name
        frame = sample_data.read_features()

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

    def _get_chunk_size(self, number_samples: int):
        max_chunk_size = 25
        chunk_size = max(1, int(number_samples / self._processing_cores))
        return min(max_chunk_size, chunk_size)

    def _read_features_table(self) -> pd.DataFrame:
        frames = []
        number_samples = len(self._sample_files_map)
        sample_data_list = list(self._sample_files_map.values())
        chunk_size = self._get_chunk_size(number_samples)
        logger.debug(f'Starting to read features table from {len(self._sample_files_map)} VCF files')
        logger.debug(f'Starting preprocessing {number_samples} samples '
                     f'with {self._processing_cores} cores and chunk size {chunk_size}')
        with mp.Pool(self._processing_cores) as pool:
            frame_vcf_masks = pool.imap_unordered(self.read_sample_data_features,
                                                  sample_data_list,
                                                  chunk_size)
            for frame_vcf_mask in frame_vcf_masks:
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

    def get_sample_files(self, sample_name: str) -> Optional[SampleData]:
        return self._sample_files_map[sample_name]

    def get_genomic_masked_region(self, sample_name: str) -> MaskedGenomicRegions:
        return self._sample_files_map[sample_name].get_mask()

    def samples_list(self) -> List[str]:
        return list(self._sample_files_map.keys())

    @classmethod
    def create(cls, sample_files_map: Dict[str, NucleotideSampleData],
               include_masked_regions: bool = True, ncores: int = 1):
        return cls(sample_files_map=sample_files_map, include_masked_regions=include_masked_regions,
                   ncores=ncores)
