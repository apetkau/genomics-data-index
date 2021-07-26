import logging
import shutil
from pathlib import Path
from typing import Tuple, Optional

import pandas as pd

from genomics_data_index.storage.MaskedGenomicRegions import MaskedGenomicRegions
from genomics_data_index.storage.io.SampleData import SampleData
from genomics_data_index.storage.io.mutation.VariationFile import VariationFile
from genomics_data_index.storage.io.mutation.VcfSnpEffAnnotationParser import VcfSnpEffAnnotationParser
from genomics_data_index.storage.model import NUCLEOTIDE_UNKNOWN_TYPE
from genomics_data_index.storage.util import TRACE_LEVEL

logger = logging.getLogger(__name__)


class NucleotideSampleData(SampleData):
    VCF_FRAME_COLUMNS = ['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'FILE', 'VARIANT_ID']

    def __init__(self, sample_name: str, vcf_file: Optional[Path], vcf_file_index: Optional[Path],
                 mask_bed_file: Optional[Path],
                 preprocessed: bool):
        super().__init__(sample_name=sample_name)
        self._vcf_file = vcf_file
        self._vcf_file_index = vcf_file_index
        self._preprocessed = preprocessed
        self._mask_bed_file = mask_bed_file
        self._snpeff_parser = VcfSnpEffAnnotationParser()

    def _preprocess_mask(self, output_dir: Path) -> Path:
        if self._mask_bed_file is None:
            mask_file = output_dir / f'{self.sample_name_persistence}.bed.gz'
            self._assert_file_not_exists(mask_file, 'Cannot preprocess data')
            mask = MaskedGenomicRegions.empty_mask()
            mask.write(mask_file)
            self._mask_bed_file = mask_file
        return self._mask_bed_file

    def _preprocess_vcf(self, output_dir: Path) -> Tuple[Path, Path]:
        new_file = output_dir / f'{self.sample_name_persistence}.vcf.gz'
        self._assert_file_not_exists(new_file, 'Cannot preprocess data')
        return VariationFile(self._vcf_file).write(new_file)

    def _do_preprocess_and_persist(self, output_dir: Path) -> SampleData:
        processed_vcf, processed_vcf_index = self._preprocess_vcf(output_dir)
        processed_mask = self._preprocess_mask(output_dir)

        return NucleotideSampleData(sample_name=self.sample_name,
                                    vcf_file=processed_vcf,
                                    vcf_file_index=processed_vcf_index,
                                    mask_bed_file=processed_mask,
                                    preprocessed=True)

    def _do_persist(self, output_dir: Path) -> SampleData:
        if self._mask_bed_file is None:
            raise Exception('mask_bed_file is None')

        new_vcf_file = output_dir / self._vcf_file.name
        new_vcf_index = output_dir / self._vcf_file_index.name
        new_mask_bed_file = output_dir / self._mask_bed_file.name

        self._assert_file_not_exists(new_vcf_file, 'Cannot persist data')
        self._assert_file_not_exists(new_vcf_index, 'Cannot persist data')
        self._assert_file_not_exists(new_mask_bed_file, 'Cannot persist data')

        logger.log(TRACE_LEVEL, f'Copying VCF and BED files to [{output_dir}] for sample [{self.sample_name}]')

        shutil.copy(self._vcf_file, new_vcf_file)
        shutil.copy(self._vcf_file_index, new_vcf_index)
        shutil.copy(self._mask_bed_file, new_mask_bed_file)

        return NucleotideSampleData(sample_name=self.sample_name,
                                    vcf_file=new_vcf_file,
                                    vcf_file_index=new_vcf_index,
                                    mask_bed_file=new_mask_bed_file,
                                    preprocessed=True)

    def is_preprocessed(self) -> bool:
        return self._preprocessed

    def _assert_file_not_exists(self, file: Path, error_prefix: str):
        if file.exists():
            raise Exception(f'{error_prefix} for sample [{self.sample_name}]: file [{file}] exists')

    def read_vcf_features(self) -> pd.DataFrame:
        vcf_file, vcf_file_index = self.get_vcf_file()
        logger.log(TRACE_LEVEL, f'Reading data for sample [{self.sample_name}] from VCF file [{vcf_file}]')
        return VariationFile(vcf_file).read_features(self.sample_name, snpeff_parser=self._snpeff_parser)

    def read_sample_data_features(self, include_masked_regions: bool = True) -> pd.DataFrame:
        vcf_features = self.read_vcf_features()
        vcf_file, vcf_file_index = self.get_vcf_file()

        if include_masked_regions:
            logger.log(TRACE_LEVEL, f'Creating unknown/missing features for sample=[{self.sample_name}]')
            frame_mask = self.get_mask().mask_to_features()
            frame_mask['SAMPLE'] = self.sample_name
            frame_mask['FILE'] = vcf_file.name
            logger.log(TRACE_LEVEL,
                       f'Combining VCF and unknown/missing (mask) dataframes for sample=[{self.sample_name}]')
            frame_vcf_mask = self.combine_vcf_mask(vcf_features, frame_mask)
        else:
            frame_vcf_mask = vcf_features

        return frame_vcf_mask

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

    def get_vcf_file(self, ignore_preprocessed=False) -> Tuple[Path, Path]:
        if ignore_preprocessed or self._preprocessed:
            return self._vcf_file, self._vcf_file_index
        else:
            raise Exception(f'VCF file for sample [{self.sample_name}] is not preprocessed: {self._vcf_file}')

    def get_mask(self) -> MaskedGenomicRegions:
        return MaskedGenomicRegions.from_file(self.get_mask_file())

    def get_mask_file(self) -> Path:
        if self._preprocessed and self._mask_bed_file is not None:
            return self._mask_bed_file
        else:
            raise Exception(f'Sample mask file is not preprocessed for sample [{self.sample_name}]')
