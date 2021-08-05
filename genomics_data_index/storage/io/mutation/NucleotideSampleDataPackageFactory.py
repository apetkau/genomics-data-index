from typing import List, Dict, Tuple
from pathlib import Path
from os import path, listdir
import logging
import pandas as pd

from genomics_data_index.storage.io.mutation.variants_processor.SerialVcfVariantsTableProcessor import \
    SerialVcfVariantsTableProcessorFactory

from genomics_data_index.storage.io.mutation.variants_processor.MultipleProcessVcfVariantsTableProcessor import \
    MultipleProcessVcfVariantsTableProcessorFactory

from genomics_data_index.storage.io.SampleDataPackageFactory import SampleDataPackageFactory
from genomics_data_index.storage.io.SampleDataPackage import SampleDataPackage
from genomics_data_index.storage.io.mutation.NucleotideSampleDataPackage import NucleotideSampleDataPackage
from genomics_data_index.storage.io.SampleFilesProcessor import SampleFilesProcessor
from genomics_data_index.storage.io.mutation.variants_processor.VcfVariantsTableProcessor import \
    VcfVariantsTableProcessorFactory
from genomics_data_index.storage.io.processor.MultipleProcessSampleFilesProcessor import \
    MultipleProcessSampleFilesProcessor
from genomics_data_index.storage.io.processor.NullSampleFilesProcessor import NullSampleFilesProcessor

logger = logging.getLogger(__file__)


class NucleotideSampleDataPackageFactory(SampleDataPackageFactory):

    def __init__(self, ncores: int, index_unknown: bool, preprocess_dir: Path):
        super().__init__()
        self._ncores = ncores
        self._preprocess_dir = preprocess_dir
        self._index_unknown = index_unknown
        self._sample_files_processor, self._variants_processor_factory = self._create_file_processors()

    def expected_input_columns(self) -> List[str]:
        return ['Sample', 'VCF', 'Mask File']

    def create_sample_vcf_mask(self, input_files_file) -> Tuple[Dict[str, Path], Dict[str, Path]]:
        logger.debug(f'Examining files listed in {input_files_file}')
        sample_vcf_map = {}
        mask_files_map = {}
        files_df = pd.read_csv(input_files_file, sep='\t')
        for index, row in files_df.iterrows():
            if row['Sample'] in sample_vcf_map:
                raise Exception(f'Error, duplicate samples {row["Sample"]} in file {input_files_file}')

            sample_vcf_map[row['Sample']] = row['VCF']
            if not pd.isna(row['Mask File']):
                mask_files_map[row['Sample']] = row['Mask File']

        logger.info(f'Found {len(sample_vcf_map)} VCFs and {len(mask_files_map)} mask files in [{input_files_file}]')

        return sample_vcf_map, mask_files_map

    def _create_file_processors(self) -> Tuple[SampleFilesProcessor, VcfVariantsTableProcessorFactory]:
        if self._ncores > 1:
            file_processor = MultipleProcessSampleFilesProcessor(preprocess_dir=Path(self._preprocess_dir),
                                                                 processing_cores=self._ncores)
            variants_processor_factory = MultipleProcessVcfVariantsTableProcessorFactory(ncores=self._ncores)
        else:
            file_processor = NullSampleFilesProcessor.instance()
            variants_processor_factory = SerialVcfVariantsTableProcessorFactory.instance()

        return file_processor, variants_processor_factory

    def create_data_package_from_snippy(self, snippy_dir: Path) -> SampleDataPackage:
        snippy_dir = Path(snippy_dir)
        sample_dirs = [snippy_dir / d for d in listdir(snippy_dir) if path.isdir(snippy_dir / d)]
        logger.debug(f'Found {len(sample_dirs)} directories in snippy_dir=[{snippy_dir}], loading files from '
                     f'these directories.')
        data_package = NucleotideSampleDataPackage.create_from_snippy(sample_dirs=sample_dirs,
                                                                      sample_files_processor=self._sample_files_processor,
                                                                      variants_processor_factory=self._variants_processor_factory,
                                                                      index_unknown_missing=self._index_unknown)
        return data_package

    def create_data_package(self, input_files_file: Path) -> SampleDataPackage:
        sample_vcf, mask_files = self.create_sample_vcf_mask(input_files_file)

        data_package = NucleotideSampleDataPackage.create_from_sequence_masks(sample_vcf_map=sample_vcf,
                                                                              masked_genomic_files_map=mask_files,
                                                                              sample_files_processor=self._sample_files_processor,
                                                                              variants_processor_factory=self._variants_processor_factory,
                                                                              index_unknown_missing=self._index_unknown)
        return data_package
