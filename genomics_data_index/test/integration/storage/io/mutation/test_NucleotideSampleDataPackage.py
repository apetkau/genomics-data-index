import re
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import cast, Dict

import pytest
from pybedtools import BedTool

from genomics_data_index.storage.io.mutation.NucleotideSampleData import NucleotideSampleData
from genomics_data_index.storage.io.mutation.NucleotideSampleDataPackage import NucleotideSampleDataPackage
from genomics_data_index.storage.io.mutation.variants_processor.MultipleProcessVcfVariantsTableProcessor import \
    MultipleProcessVcfVariantsTableProcessorFactory
from genomics_data_index.storage.io.mutation.variants_processor.SerialVcfVariantsTableProcessor import \
    SerialVcfVariantsTableProcessorFactory
from genomics_data_index.storage.io.processor.MultipleProcessSampleFilesProcessor import \
    MultipleProcessSampleFilesProcessor
from genomics_data_index.storage.io.processor.SerialSampleFilesProcessor import SerialSampleFilesProcessor
from genomics_data_index.test.integration import data_dir
from genomics_data_index.test.integration.storage.io.mutation import vcf_and_mask_files, vcf_and_bed_mask_files


def vcf_and_mixed_mask_files() -> Dict[str, Dict[str, Path]]:
    sample_vcf_map = {}
    sample_genomic_files_mask = {}

    sample_vcf_map['SampleA'] = data_dir / 'SampleA' / 'snps.vcf.gz'
    sample_genomic_files_mask['SampleA'] = data_dir / 'SampleA' / 'snps.aligned.bed.gz'
    sample_vcf_map['SampleB'] = data_dir / 'SampleB' / 'snps.vcf.gz'
    sample_genomic_files_mask['SampleB'] = data_dir / 'SampleB' / 'snps.aligned.fa'
    sample_vcf_map['SampleC'] = data_dir / 'SampleC' / 'snps.vcf.gz'
    sample_genomic_files_mask['SampleC'] = None

    return {
        'vcfs': sample_vcf_map,
        'masks': sample_genomic_files_mask
    }


def test_iter_sample_data(sample_dirs):
    with TemporaryDirectory() as tmp_file_str:
        tmp_file = Path(tmp_file_str)
        vcf_masks = vcf_and_mask_files(sample_dirs)
        file_processor = SerialSampleFilesProcessor(tmp_file)
        data_package = NucleotideSampleDataPackage.create_from_vcf_masks(sample_vcf_map=vcf_masks['vcfs'],
                                                                         masked_genomic_files_map=vcf_masks[
                                                                             'masks'],
                                                                         sample_files_processor=file_processor)

        processed_files_dict = {}
        for sample_file in data_package.iter_sample_data():
            processed_files_dict[sample_file.sample_name] = sample_file

        assert 3 == len(processed_files_dict)

        count = 0
        for sample, mask_len in [('SampleA', 437), ('SampleB', 276), ('SampleC', 329)]:
            sample_data = processed_files_dict[sample]
            assert isinstance(sample_data, NucleotideSampleData)
            sample_data = cast(NucleotideSampleData, sample_data)
            mask = sample_data.get_mask()
            assert mask_len == len(mask)
            assert {'reference'} == mask.sequence_names()
            mask_file = sample_data.get_mask_file()
            assert mask_file.parent == tmp_file
            vcf_file, vcf_index = sample_data.get_vcf_file()
            assert vcf_file.exists()
            assert vcf_index.exists()
            assert vcf_file.parent == tmp_file
            assert vcf_index.parent == tmp_file

            count += 1

        # Make sure I tested all 3 files in above loop
        assert 3 == count


def test_iter_sample_data_snippy(sample_dirs):
    with TemporaryDirectory() as tmp_file_str:
        tmp_file = Path(tmp_file_str)
        file_processor = SerialSampleFilesProcessor(tmp_file)
        expected_data = vcf_and_bed_mask_files(sample_dirs)
        data_package = NucleotideSampleDataPackage.create_from_snippy(sample_dirs,
                                                                      sample_files_processor=file_processor)

        processed_files_dict = {}
        for sample_file in data_package.iter_sample_data():
            processed_files_dict[sample_file.sample_name] = sample_file

        assert 3 == len(processed_files_dict)

        count = 0
        for sample in ['SampleA', 'SampleB', 'SampleC']:
            sample_data = processed_files_dict[sample]
            assert isinstance(sample_data, NucleotideSampleData)
            sample_data = cast(NucleotideSampleData, sample_data)
            mask = sample_data.get_mask()

            expected_mask = BedTool(str(expected_data['masks_minus_vcf'][sample]))
            assert expected_mask == mask.mask, f'expected=\n{expected_mask} != actual=\n{mask.mask}'

            assert {'reference'} == mask.sequence_names()
            mask_file = sample_data.get_mask_file()
            assert mask_file.parent == tmp_file
            vcf_file, vcf_index = sample_data.get_vcf_file()
            assert vcf_file.exists()
            assert vcf_index.exists()
            assert vcf_file.parent == tmp_file
            assert vcf_index.parent == tmp_file

            count += 1

        # Make sure I tested all 3 files in above loop
        assert 3 == count


def test_iter_sample_data_bed_masks(sample_dirs):
    with TemporaryDirectory() as tmp_file_str:
        tmp_file = Path(tmp_file_str)
        vcf_masks = vcf_and_bed_mask_files(sample_dirs)
        file_processor = SerialSampleFilesProcessor(tmp_file)
        data_package = NucleotideSampleDataPackage.create_from_vcf_masks(sample_vcf_map=vcf_masks['vcfs'],
                                                                         masked_genomic_files_map=vcf_masks[
                                                                             'masks'],
                                                                         sample_files_processor=file_processor)

        processed_files_dict = {}
        for sample_file in data_package.iter_sample_data():
            processed_files_dict[sample_file.sample_name] = sample_file

        assert 3 == len(processed_files_dict)
        assert {'SampleA', 'SampleB', 'SampleC'} == set(processed_files_dict.keys())

        count = 0
        for sample, mask_len in [('SampleA', 437), ('SampleB', 276), ('SampleC', 329)]:
            sample_data = processed_files_dict[sample]
            assert isinstance(sample_data, NucleotideSampleData)
            sample_data = cast(NucleotideSampleData, sample_data)
            mask = sample_data.get_mask()
            assert mask_len == len(mask)
            assert {'reference'} == mask.sequence_names()
            mask_file = sample_data.get_mask_file()
            assert mask_file.parent == tmp_file
            vcf_file, vcf_index = sample_data.get_vcf_file()
            assert vcf_file.exists()
            assert vcf_index.exists()
            assert vcf_file.parent == tmp_file
            assert vcf_index.parent == tmp_file

            count += 1

        # Make sure I tested all 3 files in above loop
        assert 3 == count


def test_iter_sample_data_mixed_masks(sample_dirs):
    with TemporaryDirectory() as tmp_file_str:
        tmp_file = Path(tmp_file_str)
        vcf_masks = vcf_and_mixed_mask_files()
        file_processor = SerialSampleFilesProcessor(tmp_file)
        data_package = NucleotideSampleDataPackage.create_from_vcf_masks(sample_vcf_map=vcf_masks['vcfs'],
                                                                         masked_genomic_files_map=vcf_masks[
                                                                             'masks'],
                                                                         sample_files_processor=file_processor)

        processed_files_dict = {}
        for sample_file in data_package.iter_sample_data():
            processed_files_dict[sample_file.sample_name] = sample_file

        assert 3 == len(processed_files_dict)
        assert {'SampleA', 'SampleB', 'SampleC'} == set(processed_files_dict.keys())

        count = 0
        for sample, mask_len, seq_names in [('SampleA', 437, {'reference'}), ('SampleB', 276, {'reference'}),
                                            ('SampleC', 0, set())]:
            sample_data = processed_files_dict[sample]
            assert isinstance(sample_data, NucleotideSampleData)
            sample_data = cast(NucleotideSampleData, sample_data)
            mask = sample_data.get_mask()
            assert mask_len == len(mask)
            assert seq_names == mask.sequence_names()
            mask_file = sample_data.get_mask_file()
            assert mask_file.parent == tmp_file
            vcf_file, vcf_index = sample_data.get_vcf_file()
            assert vcf_file.exists()
            assert vcf_index.exists()
            assert vcf_file.parent == tmp_file
            assert vcf_index.parent == tmp_file

            count += 1

        # Make sure I tested all 3 files in above loop
        assert 3 == count


def do_test_get_features_reader_and_features_table(data_package: NucleotideSampleDataPackage):
    # If data is not preprocessed, make sure exception is raised
    features_reader = data_package.get_features_reader()
    assert features_reader is not None
    with pytest.raises(Exception) as execinfo:
        features_reader.get_features_table()
    assert 'VCF file for sample' in str(execinfo.value)
    assert 'is not preprocessed' in str(execinfo.value)

    # Process data
    data_package_processed = data_package.process_all_data()

    # Now get features reader from new data package and verify I can read feature tables
    features_reader = data_package_processed.get_features_reader()
    features_df = features_reader.get_features_table()
    assert 1170 == len(features_df)
    assert {'SampleA', 'SampleB', 'SampleC'} == set(features_df['SAMPLE'].tolist())


def test_get_features_reader_and_features_table_serial_variants_processor(sample_dirs):
    with TemporaryDirectory() as tmp_file_str:
        tmp_file = Path(tmp_file_str)
        vcf_masks = vcf_and_mask_files(sample_dirs)
        file_processor = SerialSampleFilesProcessor(tmp_file)
        variants_processor_factory = SerialVcfVariantsTableProcessorFactory.instance()
        data_package = NucleotideSampleDataPackage.create_from_vcf_masks(sample_vcf_map=vcf_masks['vcfs'],
                                                                         masked_genomic_files_map=vcf_masks[
                                                                             'masks'],
                                                                         variants_processor_factory=variants_processor_factory,
                                                                         sample_files_processor=file_processor)
        do_test_get_features_reader_and_features_table(data_package)


def test_get_features_reader_and_features_table_parallel_variants_processor(sample_dirs):
    with TemporaryDirectory() as tmp_file_str:
        tmp_file = Path(tmp_file_str)
        vcf_masks = vcf_and_mask_files(sample_dirs)
        file_processor = SerialSampleFilesProcessor(tmp_file)
        variants_processor_factory = MultipleProcessVcfVariantsTableProcessorFactory(ncores=2)
        data_package = NucleotideSampleDataPackage.create_from_vcf_masks(sample_vcf_map=vcf_masks['vcfs'],
                                                                         masked_genomic_files_map=vcf_masks[
                                                                             'masks'],
                                                                         variants_processor_factory=variants_processor_factory,
                                                                         sample_files_processor=file_processor)
        do_test_get_features_reader_and_features_table(data_package)


def test_persisted_sample_data_file_names(sample_dirs):
    def split_all_ext(file: Path) -> str:
        return file.name.split('.')[0]

    with TemporaryDirectory() as tmp_file_str:
        tmp_file = Path(tmp_file_str)
        vcf_masks = vcf_and_mask_files(sample_dirs)
        file_processor = SerialSampleFilesProcessor(tmp_file)
        data_package = NucleotideSampleDataPackage.create_from_vcf_masks(sample_vcf_map=vcf_masks['vcfs'],
                                                                         masked_genomic_files_map=vcf_masks[
                                                                             'masks'],
                                                                         sample_files_processor=file_processor)

        processed_data_package = cast(NucleotideSampleDataPackage, data_package.process_all_data())
        processed_files_dict = processed_data_package.get_sample_data()
        sample_data = cast(NucleotideSampleData, processed_files_dict['SampleA'])
        mask_file = sample_data.get_mask_file()
        vcf_file, vcf_index = sample_data.get_vcf_file()

        assert sample_data.sample_name_persistence == split_all_ext(mask_file)
        assert sample_data.sample_name_persistence == split_all_ext(vcf_file)
        assert sample_data.sample_name_persistence == split_all_ext(vcf_index)
        assert 32 == len(sample_data.sample_name_persistence)
        assert bool(re.match(r'^[0-9a-f]*$', sample_data.sample_name_persistence))


def test_with_serial_sample_files_processor(sample_dirs):
    with TemporaryDirectory() as tmp_file_str:
        tmp_file = Path(tmp_file_str)
        vcf_masks = vcf_and_mask_files(sample_dirs)
        file_processor = SerialSampleFilesProcessor(tmp_file)
        data_package = NucleotideSampleDataPackage.create_from_vcf_masks(sample_vcf_map=vcf_masks['vcfs'],
                                                                         masked_genomic_files_map=vcf_masks[
                                                                             'masks'],
                                                                         sample_files_processor=file_processor)

        processed_data_package = cast(NucleotideSampleDataPackage, data_package.process_all_data())
        processed_files_dict = processed_data_package.get_sample_data()

        assert 3 == len(processed_files_dict)

        sample_data = processed_files_dict['SampleA']
        assert isinstance(sample_data, NucleotideSampleData)
        sample_data = cast(NucleotideSampleData, sample_data)
        mask = sample_data.get_mask()
        assert 437 == len(mask)
        assert {'reference'} == mask.sequence_names()
        assert sample_data.get_mask_file().parent == tmp_file
        vcf_file, vcf_index = sample_data.get_vcf_file()
        assert vcf_file.exists()
        assert vcf_index.exists()
        assert vcf_file.parent == tmp_file
        assert vcf_index.parent == tmp_file

        sample_data = processed_files_dict['SampleB']
        assert isinstance(sample_data, NucleotideSampleData)
        sample_data = cast(NucleotideSampleData, sample_data)
        mask = sample_data.get_mask()
        assert 276 == len(mask)
        assert {'reference'} == mask.sequence_names()
        assert sample_data.get_mask_file().parent == tmp_file
        vcf_file, vcf_index = sample_data.get_vcf_file()
        assert vcf_file.exists()
        assert vcf_index.exists()
        assert vcf_file.parent == tmp_file
        assert vcf_index.parent == tmp_file

        sample_data = processed_files_dict['SampleC']
        assert isinstance(sample_data, NucleotideSampleData)
        sample_data = cast(NucleotideSampleData, sample_data)
        mask = sample_data.get_mask()
        assert 329 == len(mask)
        assert {'reference'} == mask.sequence_names()
        assert sample_data.get_mask_file().parent == tmp_file
        vcf_file, vcf_index = sample_data.get_vcf_file()
        assert vcf_file.exists()
        assert vcf_index.exists()
        assert vcf_file.parent == tmp_file
        assert vcf_index.parent == tmp_file


def test_with_multiprocess_sample_files_processor(sample_dirs):
    with TemporaryDirectory() as tmp_file_str:
        tmp_file = Path(tmp_file_str)
        vcf_masks = vcf_and_mask_files(sample_dirs)
        file_processor = MultipleProcessSampleFilesProcessor(tmp_file, processing_cores=2)
        data_package = NucleotideSampleDataPackage.create_from_vcf_masks(sample_vcf_map=vcf_masks['vcfs'],
                                                                         masked_genomic_files_map=vcf_masks[
                                                                             'masks'],
                                                                         sample_files_processor=file_processor)

        processed_data_package = cast(NucleotideSampleDataPackage, data_package.process_all_data())
        processed_files_dict = processed_data_package.get_sample_data()

        assert 3 == len(processed_files_dict)

        sample_data = processed_files_dict['SampleA']
        assert isinstance(sample_data, NucleotideSampleData)
        sample_data = cast(NucleotideSampleData, sample_data)
        mask = sample_data.get_mask()
        assert 437 == len(mask)
        assert {'reference'} == mask.sequence_names()
        assert sample_data.get_mask_file().parent == tmp_file
        vcf_file, vcf_index = sample_data.get_vcf_file()
        assert vcf_file.exists()
        assert vcf_index.exists()
        assert vcf_file.parent == tmp_file
        assert vcf_index.parent == tmp_file

        sample_data = processed_files_dict['SampleB']
        assert isinstance(sample_data, NucleotideSampleData)
        sample_data = cast(NucleotideSampleData, sample_data)
        mask = sample_data.get_mask()
        assert 276 == len(mask)
        assert {'reference'} == mask.sequence_names()
        assert sample_data.get_mask_file().parent == tmp_file
        vcf_file, vcf_index = sample_data.get_vcf_file()
        assert vcf_file.exists()
        assert vcf_index.exists()
        assert vcf_file.parent == tmp_file
        assert vcf_index.parent == tmp_file

        sample_data = processed_files_dict['SampleC']
        assert isinstance(sample_data, NucleotideSampleData)
        sample_data = cast(NucleotideSampleData, sample_data)
        mask = sample_data.get_mask()
        assert 329 == len(mask)
        assert {'reference'} == mask.sequence_names()
        assert sample_data.get_mask_file().parent == tmp_file
        vcf_file, vcf_index = sample_data.get_vcf_file()
        assert vcf_file.exists()
        assert vcf_index.exists()
        assert vcf_file.parent == tmp_file
        assert vcf_index.parent == tmp_file


def test_with_null_sample_files_processor(sample_dirs):
    vcf_masks = vcf_and_mask_files(sample_dirs)
    data_package = NucleotideSampleDataPackage.create_from_vcf_masks(sample_vcf_map=vcf_masks['vcfs'],
                                                                     masked_genomic_files_map=vcf_masks[
                                                                         'masks'])

    processed_data_package = cast(NucleotideSampleDataPackage, data_package.process_all_data())
    processed_files_dict = processed_data_package.get_sample_data()

    assert 3 == len(processed_files_dict)
    sample_data = processed_files_dict['SampleA']
    assert isinstance(sample_data, NucleotideSampleData)
    sample_data = cast(NucleotideSampleData, sample_data)

    with pytest.raises(Exception) as execinfo:
        sample_data.get_mask()
    assert 'Sample mask file is not preprocessed for sample' in str(execinfo.value)

    with pytest.raises(Exception) as execinfo:
        sample_data.get_vcf_file()
    assert 'VCF file for sample' in str(execinfo.value)
