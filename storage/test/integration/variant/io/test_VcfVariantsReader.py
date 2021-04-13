from os import path, listdir
from pathlib import Path
from typing import List

import pytest

from storage.test.integration.variant import data_dir
from storage.test.integration.variant import data_dir_empty
from storage.variant.io.mutation.VcfVariantsReader import VcfVariantsReader


@pytest.fixture
def sample_dirs() -> List[Path]:
    return [data_dir / d for d in listdir(data_dir) if path.isdir(data_dir / d)]


@pytest.fixture
def sample_dirs_empty() -> List[Path]:
    return [data_dir_empty / d for d in listdir(data_dir_empty) if path.isdir(data_dir_empty / d)]


def variants_reader_internal(sample_dirs):
    sample_vcf_map = {}
    sample_genomic_files_mask = {}
    for d in sample_dirs:
        sample_name = path.basename(d)
        vcf_file = Path(d, 'snps.vcf.gz')
        genomic_file_mask = Path(d, 'snps.aligned.fa')

        sample_vcf_map[sample_name] = vcf_file
        sample_genomic_files_mask[sample_name] = genomic_file_mask

    return VcfVariantsReader(sample_vcf_map=sample_vcf_map,
                             masked_genomic_files_map=sample_genomic_files_mask)


@pytest.fixture
def variants_reader(sample_dirs) -> VcfVariantsReader:
    return variants_reader_internal(sample_dirs)


@pytest.fixture
def variants_reader_empty(sample_dirs_empty) -> VcfVariantsReader:
    return variants_reader_internal(sample_dirs_empty)


@pytest.fixture
def variants_reader_empty_masks(sample_dirs) -> VcfVariantsReader:
    sample_vcf_map = {}
    for d in sample_dirs:
        sample_name = path.basename(d)
        vcf_file = Path(d, 'snps.vcf.gz')

        sample_vcf_map[sample_name] = vcf_file

    return VcfVariantsReader(sample_vcf_map=sample_vcf_map)


def test_get_variants_table(variants_reader):
    df = variants_reader.get_features_table()

    assert 129 == len(df), 'Data has incorrect length'
    assert {'SampleA', 'SampleB', 'SampleC'} == set(df['SAMPLE'].tolist()), 'Incorrect sample names'

    assert ['snp'] == df.loc[(df['SAMPLE'] == 'SampleA') & (df['POS'] == 293), 'TYPE'].tolist()
    assert ['del'] == df.loc[(df['SAMPLE'] == 'SampleA') & (df['POS'] == 302), 'TYPE'].tolist()
    assert ['ins'] == df.loc[(df['SAMPLE'] == 'SampleA') & (df['POS'] == 324), 'TYPE'].tolist()
    assert ['ins'] == df.loc[(df['SAMPLE'] == 'SampleA') & (df['POS'] == 374), 'TYPE'].tolist()
    assert ['complex'] == df.loc[(df['SAMPLE'] == 'SampleA') & (df['POS'] == 461), 'TYPE'].tolist()
    assert ['complex'] == df.loc[(df['SAMPLE'] == 'SampleB') & (df['POS'] == 1325), 'TYPE'].tolist()
    assert ['complex'] == df.loc[(df['SAMPLE'] == 'SampleC') & (df['POS'] == 1984), 'TYPE'].tolist()


def test_get_genomic_masks(variants_reader):
    genomic_masks = variants_reader.get_genomic_masked_regions()

    assert {'SampleA', 'SampleB', 'SampleC'} == set(genomic_masks.keys()), 'Incorrect samples'
    assert 437 == len(genomic_masks['SampleA'])
    assert 276 == len(genomic_masks['SampleB'])
    assert 329 == len(genomic_masks['SampleC'])

    assert {'reference'} == genomic_masks['SampleA'].sequence_names()
    assert {'reference'} == genomic_masks['SampleB'].sequence_names()
    assert {'reference'} == genomic_masks['SampleC'].sequence_names()


def test_get_genomic_masks_empty(variants_reader_empty_masks):
    genomic_masks = variants_reader_empty_masks.get_genomic_masked_regions()
    assert {'SampleA', 'SampleB', 'SampleC'} == set(genomic_masks.keys())
    assert genomic_masks['SampleA'].is_empty()


def test_get_samples_list(variants_reader):
    assert {'SampleA', 'SampleB', 'SampleC'} == set(variants_reader.samples_list())


def test_get_variants_table_empty(variants_reader_empty):
    df = variants_reader_empty.get_features_table()

    assert 0 == len(df), 'Data has incorrect length'


def test_get_or_create_feature_file(variants_reader):
    file = variants_reader.get_or_create_feature_file('SampleA')
    assert file.exists()
    assert 'SampleA' in str(file)
