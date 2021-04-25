from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Dict

import pytest

from genomics_data_index.storage.io.mutation.NucleotideSampleDataSequenceMask import NucleotideSampleDataSequenceMask
from genomics_data_index.test.integration import regular_vcf_dir, data_dir
from genomics_data_index.test.integration.storage.io import read_vcf_df


@pytest.fixture
def files_map():
    return {
        'SampleA': {
            'vcf': Path(regular_vcf_dir, 'SampleA.vcf.gz'),
            'mask': Path(data_dir, 'SampleA', 'snps.aligned.fa'),
        },
        'SampleB': {
            'vcf': Path(regular_vcf_dir, 'SampleB.vcf.gz'),
            'mask': Path(data_dir, 'SampleB', 'snps.aligned.fa'),
        },
        'SampleC': {
            'vcf': Path(regular_vcf_dir, 'SampleC.vcf.gz'),
            'mask': Path(data_dir, 'SampleC', 'snps.aligned.fa'),
        }
    }


def test_get_vcf_file_no_preprocess(files_map: Dict[str, Dict[str, Path]]):
    sample_files = NucleotideSampleDataSequenceMask.create(sample_name='SampleA',
                                                           vcf_file=files_map['SampleA']['vcf'],
                                                           sample_mask_sequence=files_map['SampleA']['mask'])
    vcf_file, index_file = sample_files.get_vcf_file(ignore_preprocessed=True)
    assert vcf_file.exists()
    assert index_file is None
    df = read_vcf_df(vcf_file)
    assert len(df[df['TYPE'].isna()]) > 0


def test_get_sample_mask_no_preprocess(files_map: Dict[str, Dict[str, Path]]):
    sample_files = NucleotideSampleDataSequenceMask.create(sample_name='SampleA',
                                                           vcf_file=files_map['SampleA']['vcf'],
                                                           sample_mask_sequence=files_map['SampleA']['mask'])
    with pytest.raises(Exception) as execinfo:
        sample_files.get_mask()
    assert 'Sample mask file is not preprocessed for sample' in str(execinfo.value)


def test_get_empty_sample_mask_no_preprocess(files_map: Dict[str, Dict[str, Path]]):
    sample_files = NucleotideSampleDataSequenceMask.create(sample_name='SampleA',
                                                           vcf_file=files_map['SampleA']['vcf'],
                                                           sample_mask_sequence=None)
    with pytest.raises(Exception) as execinfo:
        sample_files.get_mask()
    assert 'Sample mask file is not preprocessed for sample' in str(execinfo.value)


def test_get_sample_mask_after_preprocess(files_map: Dict[str, Dict[str, Path]]):
    with TemporaryDirectory() as tmp_dir_str:
        tmp_dir = Path(tmp_dir_str)
        sample_files = NucleotideSampleDataSequenceMask.create(sample_name='SampleA',
                                                               vcf_file=files_map['SampleA']['vcf'],
                                                               sample_mask_sequence=files_map['SampleA']['mask'])
        assert not sample_files.is_preprocessed()
        processed_sample_files = sample_files.persist(tmp_dir)
        assert processed_sample_files.is_preprocessed()
        mask = processed_sample_files.get_mask()
        assert 437 == len(mask)
        assert {'reference'} == mask.sequence_names()
        mask_file = processed_sample_files.get_mask_file()
        assert mask_file.exists()
        assert mask_file.parent == tmp_dir


def test_get_sample_empty_mask_after_preprocess(files_map: Dict[str, Dict[str, Path]]):
    with TemporaryDirectory() as tmp_dir_str:
        tmp_dir = Path(tmp_dir_str)
        sample_files = NucleotideSampleDataSequenceMask.create(sample_name='SampleA',
                                                               vcf_file=files_map['SampleA']['vcf'])
        assert not sample_files.is_preprocessed()
        processed_sample_files = sample_files.persist(tmp_dir)
        assert processed_sample_files.is_preprocessed()
        mask_file = processed_sample_files.get_mask_file()
        assert mask_file.exists()
        assert mask_file.parent == tmp_dir
        mask = processed_sample_files.get_mask()
        assert mask.is_empty()


def test_get_sample_mask_double_preprocess(files_map: Dict[str, Dict[str, Path]]):
    with TemporaryDirectory() as tmp_dir_str:
        tmp_dir = Path(tmp_dir_str)
        sample_files = NucleotideSampleDataSequenceMask.create(sample_name='SampleA',
                                                               vcf_file=files_map['SampleA']['vcf'],
                                                               sample_mask_sequence=files_map['SampleA']['mask'])
        assert not sample_files.is_preprocessed()
        processed_sample_files = sample_files.persist(tmp_dir)
        assert processed_sample_files.is_preprocessed()
        mask = processed_sample_files.get_mask()
        assert 437 == len(mask)
        assert {'reference'} == mask.sequence_names()
        mask_file = processed_sample_files.get_mask_file()

        with TemporaryDirectory() as tmp_dir_2_str:
            tmp_dir_2 = Path(tmp_dir_2_str)
            processed_sample_files_2 = processed_sample_files.persist(tmp_dir_2)
            assert processed_sample_files_2.is_preprocessed()
            mask_2 = processed_sample_files_2.get_mask()
            assert 437 == len(mask)
            assert {'reference'} == mask_2.sequence_names()
            mask_file_2 = processed_sample_files_2.get_mask_file()

            assert mask_file.exists()
            assert mask_file_2.exists()
            assert mask_file.parent == tmp_dir
            assert mask_file_2.parent == tmp_dir_2


def test_get_vcf_file_no_preprocess_exception(files_map: Dict[str, Dict[str, Path]]):
    sample_files = NucleotideSampleDataSequenceMask.create(sample_name='SampleA',
                                                           vcf_file=files_map['SampleA']['vcf'],
                                                           sample_mask_sequence=files_map['SampleA']['mask'])
    with pytest.raises(Exception) as execinfo:
        sample_files.get_vcf_file()
    assert 'VCF file for sample' in str(execinfo.value)


def test_get_vcf_file_with_preprocess(files_map: Dict[str, Dict[str, Path]]):
    with TemporaryDirectory() as tmp_dir_str:
        tmp_dir = Path(tmp_dir_str)
        sample_files = NucleotideSampleDataSequenceMask.create(sample_name='SampleA',
                                                               vcf_file=files_map['SampleA']['vcf'],
                                                               sample_mask_sequence=files_map['SampleA']['mask'])
        assert not sample_files.is_preprocessed()
        with pytest.raises(Exception) as execinfo:
            sample_files.get_vcf_file()
        assert 'VCF file for sample' in str(execinfo.value)

        processed_sample_files = sample_files.persist(tmp_dir)
        assert processed_sample_files.is_preprocessed()

        vcf_file, index_file = processed_sample_files.get_vcf_file()
        assert vcf_file.exists()
        assert index_file.exists()
        assert vcf_file.parent == tmp_dir
        assert index_file.parent == tmp_dir

        df = read_vcf_df(vcf_file)
        assert len(df[df['TYPE'].isna()]) == 0
        assert 'SNP' == df[df['POS'] == 293]['TYPE'].tolist()[0]
        assert 'INDEL' == df[df['POS'] == 302]['TYPE'].tolist()[0]
        assert 'INDEL' == df[df['POS'] == 324]['TYPE'].tolist()[0]
        assert 'INDEL' == df[df['POS'] == 374]['TYPE'].tolist()[0]
        assert 'OTHER' == df[df['POS'] == 461]['TYPE'].tolist()[0]
        assert 'SNP' == df[df['POS'] == 506]['TYPE'].tolist()[0]


def test_get_vcf_file_double_preprocess(files_map: Dict[str, Dict[str, Path]]):
    with TemporaryDirectory() as tmp_dir_str:
        tmp_dir = Path(tmp_dir_str)
        sample_files = NucleotideSampleDataSequenceMask.create(sample_name='SampleA',
                                                               vcf_file=files_map['SampleA']['vcf'],
                                                               sample_mask_sequence=files_map['SampleA']['mask'])
        assert not sample_files.is_preprocessed()
        processed_sample_files = sample_files.persist(tmp_dir)
        assert processed_sample_files.is_preprocessed()
        vcf_file, index_file = processed_sample_files.get_vcf_file()
        assert vcf_file.exists()
        assert index_file.exists()

        with TemporaryDirectory() as tmp_dir_2_str:
            tmp_dir_2 = Path(tmp_dir_2_str)
            processed_sample_files_2 = processed_sample_files.persist(tmp_dir_2)
            assert processed_sample_files_2.is_preprocessed()
            vcf_file_2, index_file_2 = processed_sample_files_2.get_vcf_file()
            assert vcf_file_2.exists()
            assert index_file_2.exists()

            assert vcf_file.parent == tmp_dir
            assert index_file.parent == tmp_dir
            assert vcf_file_2.parent == tmp_dir_2
            assert index_file_2.parent == tmp_dir_2
