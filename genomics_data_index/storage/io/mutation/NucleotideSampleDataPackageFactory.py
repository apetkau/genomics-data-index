import abc
import math
import logging
from os import path, listdir
from pathlib import Path
from typing import List, Dict, Tuple, Generator

import pandas as pd

from genomics_data_index.storage.io.SampleDataPackage import SampleDataPackage
from genomics_data_index.storage.io.SampleDataPackageFactory import SampleDataPackageFactory
from genomics_data_index.storage.io.SampleFilesProcessor import SampleFilesProcessor
from genomics_data_index.storage.io.mutation.NucleotideSampleDataPackage import NucleotideSampleDataPackage
from genomics_data_index.storage.io.mutation.variants_processor.MultipleProcessVcfVariantsTableProcessor import \
    MultipleProcessVcfVariantsTableProcessorFactory
from genomics_data_index.storage.io.mutation.variants_processor.SerialVcfVariantsTableProcessor import \
    SerialVcfVariantsTableProcessorFactory
from genomics_data_index.storage.io.mutation.variants_processor.VcfVariantsTableProcessor import \
    VcfVariantsTableProcessorFactory
from genomics_data_index.storage.io.processor.MultipleProcessSampleFilesProcessor import \
    MultipleProcessSampleFilesProcessor
from genomics_data_index.storage.io.processor.SerialSampleFilesProcessor import SerialSampleFilesProcessor
from genomics_data_index.storage.util.ListSliceIter import ListSliceIter
from genomics_data_index.storage.util.SamplesProgressLogger import SamplesProgressLogger

logger = logging.getLogger(__file__)


class NucleotideSampleDataPackageFactory(SampleDataPackageFactory, abc.ABC):

    def __init__(self, ncores: int, index_unknown: bool, preprocess_dir: Path):
        super().__init__()
        self._ncores = ncores
        self._preprocess_dir = preprocess_dir
        self._index_unknown = index_unknown
        self._sample_files_processor, self._variants_processor_factory = self._create_file_processors()

    @abc.abstractmethod
    def number_samples(self) -> int:
        pass

    def number_batches(self, batch_size: int) -> int:
        return (int)(math.ceil(self.number_samples()/batch_size))

    def _create_file_processors(self) -> Tuple[SampleFilesProcessor, VcfVariantsTableProcessorFactory]:
        if self._ncores > 1:
            file_processor = MultipleProcessSampleFilesProcessor(preprocess_dir=Path(self._preprocess_dir),
                                                                 processing_cores=self._ncores)
            variants_processor_factory = MultipleProcessVcfVariantsTableProcessorFactory(ncores=self._ncores)
        else:
            file_processor = SerialSampleFilesProcessor(preprocess_dir=self._preprocess_dir)
            variants_processor_factory = SerialVcfVariantsTableProcessorFactory.instance()

        return file_processor, variants_processor_factory


class NucleotideInputFilesSampleDataPackageFactory(NucleotideSampleDataPackageFactory):

    def __init__(self, ncores: int, index_unknown: bool, preprocess_dir: Path,
                 input_files_file: Path):
        super().__init__(ncores=ncores, index_unknown=index_unknown, preprocess_dir=preprocess_dir)
        self._input_files_file = input_files_file
        self._sample_vcf, self._mask_files = self.create_sample_vcf_mask(input_files_file)

    def number_samples(self) -> int:
        return len(self._sample_vcf)

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

    def create_data_package_iter(self, batch_size: int = 100) -> Generator[SampleDataPackage, None, None]:
        sample_names = list(self._sample_vcf.keys())
        number_batches = self.number_batches(batch_size)
        SamplesProgressLogger.set_total_batches(number_batches)
        sample_names_iter = ListSliceIter(sample_names, batch_size)
        current_batch = 1
        for sample_names_slice in sample_names_iter.islice():
            sample_vcfs_slice = {s: self._sample_vcf[s] for s in sample_names_slice}
            mask_files_slice = {s: self._mask_files[s] for s in sample_names_slice if s in self._mask_files}
            SamplesProgressLogger.set_current_batch(current_batch)
            logger.info(f'Starting batch {current_batch}/{number_batches}')
            yield self._create_data_package(sample_vcfs=sample_vcfs_slice, sample_mask_files=mask_files_slice)
            logger.info(f'Finished batch {current_batch}/{number_batches}')
            current_batch = current_batch + 1
        SamplesProgressLogger.unset_batch_mode()
        logger.debug(f'Finished creating data packages for all {self.number_samples()} samples')

    def _create_data_package(self, sample_vcfs: Dict[str, Path],
                             sample_mask_files: Dict[str, Path]) -> SampleDataPackage:
        return NucleotideSampleDataPackage.create_from_sequence_masks(sample_vcf_map=sample_vcfs,
                                                                      masked_genomic_files_map=sample_mask_files,
                                                                      sample_files_processor=self._sample_files_processor,
                                                                      variants_processor_factory=self._variants_processor_factory,
                                                                      index_unknown_missing=self._index_unknown)

    def create_data_package(self) -> SampleDataPackage:
        sample_vcf, mask_files = self.create_sample_vcf_mask(self._input_files_file)
        return self._create_data_package(sample_vcfs=sample_vcf, sample_mask_files=mask_files)


class NucleotideSnippySampleDataPackageFactory(NucleotideSampleDataPackageFactory):

    def __init__(self, ncores: int, index_unknown: bool, preprocess_dir: Path,
                 snippy_dir: Path):
        super().__init__(ncores=ncores, index_unknown=index_unknown, preprocess_dir=preprocess_dir)
        self._snippy_dir = snippy_dir
        self._sample_dirs = [snippy_dir / d for d in listdir(snippy_dir) if path.isdir(snippy_dir / d)]

    def number_samples(self) -> int:
        return len(self._sample_dirs)

    def create_data_package_iter(self, batch_size: int = 100) -> Generator[SampleDataPackage, None, None]:
        number_batches = self.number_batches(batch_size)
        SamplesProgressLogger.set_total_batches(number_batches)
        sample_dirs_iter = ListSliceIter(self._sample_dirs, batch_size)
        current_batch = 1
        for sample_dirs in sample_dirs_iter.islice():
            SamplesProgressLogger.set_current_batch(current_batch)
            logger.info(f'Starting batch {current_batch}/{number_batches}')
            yield self._create_data_package(sample_dirs=sample_dirs)
            logger.info(f'Finished batch {current_batch}/{number_batches}')
            current_batch = current_batch + 1
        SamplesProgressLogger.unset_batch_mode()
        logger.debug(f'Finished creating data packages for all {self.number_samples()} samples')

    def _create_data_package(self, sample_dirs: List[Path]) -> SampleDataPackage:
        data_package = NucleotideSampleDataPackage.create_from_snippy(sample_dirs=sample_dirs,
                                                                      sample_files_processor=self._sample_files_processor,
                                                                      variants_processor_factory=self._variants_processor_factory,
                                                                      index_unknown_missing=self._index_unknown)
        return data_package

    def create_data_package(self) -> SampleDataPackage:
        logger.debug(
            f'Found {len(self._sample_dirs)} directories in snippy_dir=[{self._snippy_dir}], loading files from '
            f'these directories.')
        return self._create_data_package(sample_dirs=self._sample_dirs)
