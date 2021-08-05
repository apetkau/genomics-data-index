from typing import List, Dict, Tuple, Generator
import abc
from pathlib import Path
from os import path, listdir
import logging
import pandas as pd
from genomics_data_index.storage.util.SamplesProgressLogger import SamplesProgressLogger

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
from genomics_data_index.storage.io.processor.SerialSampleFilesProcessor import SerialSampleFilesProcessor
from genomics_data_index.storage.util.ListSliceIter import ListSliceIter

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

    def _create_file_processors(self) -> Tuple[SampleFilesProcessor, VcfVariantsTableProcessorFactory]:
        if self._ncores > 1:
            file_processor = MultipleProcessSampleFilesProcessor(preprocess_dir=Path(self._preprocess_dir),
                                                                 processing_cores=self._ncores)
            variants_processor_factory = MultipleProcessVcfVariantsTableProcessorFactory(ncores=self._ncores)
        else:
            file_processor = SerialSampleFilesProcessor(preprocess_dir=self._preprocess_dir)
            variants_processor_factory = SerialVcfVariantsTableProcessorFactory.instance()

        return file_processor, variants_processor_factory

    def create_data_package_iter(self, batch_size: int = 100) -> Generator[SampleDataPackage, None, None]:
        pass


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

    def create_data_package(self) -> SampleDataPackage:
        sample_vcf, mask_files = self.create_sample_vcf_mask(self._input_files_file)

        data_package = NucleotideSampleDataPackage.create_from_sequence_masks(sample_vcf_map=sample_vcf,
                                                                              masked_genomic_files_map=mask_files,
                                                                              sample_files_processor=self._sample_files_processor,
                                                                              variants_processor_factory=self._variants_processor_factory,
                                                                              index_unknown_missing=self._index_unknown)
        return data_package


class NucleotideSnippySampleDataPackageFactory(NucleotideSampleDataPackageFactory):

    def __init__(self, ncores: int, index_unknown: bool, preprocess_dir: Path,
                 snippy_dir: Path):
        super().__init__(ncores=ncores, index_unknown=index_unknown, preprocess_dir=preprocess_dir)
        self._snippy_dir = snippy_dir
        self._sample_dirs = [snippy_dir / d for d in listdir(snippy_dir) if path.isdir(snippy_dir / d)]

    def number_samples(self) -> int:
        return len(self._sample_dirs)

    def create_data_package_iter(self, batch_size: int = 100) -> Generator[SampleDataPackage, None, None]:
        number_samples = self.number_samples()
        progress_logger = SamplesProgressLogger(stage_name='Read Samples', stage_number=0, total_samples=number_samples)
        progress_logger.update_progress(0)
        sample_dirs_iter = ListSliceIter(self._sample_dirs, batch_size)
        samples_processed = 0
        for sample_dirs in sample_dirs_iter.islice():
            yield self._create_data_package(sample_dirs=sample_dirs)
            samples_processed = samples_processed + len(sample_dirs)
            progress_logger.update_progress(samples_processed)
        logger.debug(f'Finished creating data packages for all {number_samples} samples')

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
