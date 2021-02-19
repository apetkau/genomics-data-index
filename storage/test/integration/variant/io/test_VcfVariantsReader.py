from typing import List
import pytest
from os import path, listdir
from pathlib import Path

from storage.variant.io.VcfVariantsReader import VcfVariantsReader
from storage.variant.CoreBitMask import CoreBitMask
from storage.test.integration.variant import data_dir
from storage.test.integration.variant import data_dir_empty


@pytest.fixture
def sample_dirs() -> List[Path]:
    return [data_dir / d for d in listdir(data_dir) if path.isdir(data_dir / d)]


@pytest.fixture
def sample_dirs_empty() -> List[Path]:
    return [data_dir_empty / d for d in listdir(data_dir_empty) if path.isdir(data_dir_empty / d)]


def variants_reader_internal(sample_dirs):
    sample_vcf_map = {}
    sample_core_masks_map = {}
    for d in sample_dirs:
        sample_name = path.basename(d)
        vcf_file = Path(d, 'snps.vcf.gz')
        core_mask_file = Path(d, 'snps.aligned.fa')

        sample_vcf_map[sample_name] = vcf_file
        sample_core_masks_map[sample_name] = core_mask_file

    return VcfVariantsReader(sample_vcf_map=sample_vcf_map,
                             core_mask_files_map=sample_core_masks_map)


@pytest.fixture
def variants_reader(sample_dirs) -> VcfVariantsReader:
    return variants_reader_internal(sample_dirs)


@pytest.fixture
def variants_reader_empty(sample_dirs_empty) -> VcfVariantsReader:
    return variants_reader_internal(sample_dirs_empty)


@pytest.fixture
def variants_reader_empty_core_masks(sample_dirs) -> VcfVariantsReader:
    sample_vcf_map = {}
    for d in sample_dirs:
        sample_name = path.basename(d)
        vcf_file = Path(d, 'snps.vcf.gz')

        sample_vcf_map[sample_name] = vcf_file

    empty_mask = {
        'reference': CoreBitMask.empty_mask(5180)
    }

    return VcfVariantsReader(sample_vcf_map=sample_vcf_map,
                             empty_core_mask=empty_mask)


def test_get_variants_table(variants_reader):
    df = variants_reader.get_variants_table()

    assert 129 == len(df), 'Data has incorrect length'
    assert {'SampleA', 'SampleB', 'SampleC'} == set(df['SAMPLE'].tolist()), 'Incorrect sample names'

    assert ['snp'] == df.loc[(df['SAMPLE'] == 'SampleA') & (df['POS'] == 293), 'TYPE'].tolist()
    assert ['del'] == df.loc[(df['SAMPLE'] == 'SampleA') & (df['POS'] == 302), 'TYPE'].tolist()
    assert ['ins'] == df.loc[(df['SAMPLE'] == 'SampleA') & (df['POS'] == 324), 'TYPE'].tolist()
    assert ['ins'] == df.loc[(df['SAMPLE'] == 'SampleA') & (df['POS'] == 374), 'TYPE'].tolist()
    assert ['complex'] == df.loc[(df['SAMPLE'] == 'SampleA') & (df['POS'] == 461), 'TYPE'].tolist()
    assert ['complex'] == df.loc[(df['SAMPLE'] == 'SampleB') & (df['POS'] == 1325), 'TYPE'].tolist()
    assert ['complex'] == df.loc[(df['SAMPLE'] == 'SampleC') & (df['POS'] == 1984), 'TYPE'].tolist()


def test_get_core_masks(variants_reader):
    core_masks = variants_reader.get_core_masks()

    assert {'SampleA', 'SampleB', 'SampleC'} == set(core_masks.keys()), 'Incorrect samples'
    assert 4743 == core_masks['SampleA']['reference'].core_length(), 'Incorrect core length for SampleA'
    assert 4904 == core_masks['SampleB']['reference'].core_length(), 'Incorrect core length for SampleB'
    assert 4851 == core_masks['SampleC']['reference'].core_length(), 'Incorrect core length for SampleC'


def test_get_core_masks_empty(variants_reader_empty_core_masks):
    core_masks = variants_reader_empty_core_masks.get_core_masks()
    assert {'SampleA', 'SampleB', 'SampleC'} == set(core_masks.keys())
    assert core_masks['SampleA']['reference'].is_empty()
    assert 5180 == core_masks['SampleA']['reference'].core_length()


def test_get_samples_list(variants_reader):
    assert {'SampleA', 'SampleB', 'SampleC'} == set(variants_reader.samples_list())


def test_get_variants_table_empty(variants_reader_empty):
    df = variants_reader_empty.get_variants_table()

    assert 0 == len(df), 'Data has incorrect length'
