from os import path, listdir
from typing import Dict, Any

import pytest

from storage.test.integration.variant import data_dir
from storage.test.integration.variant import data_dir_empty
from storage.variant.io.mutation.SnippyVariantsReader import SnippyVariantsReader

sample_dirs = [data_dir / d for d in listdir(data_dir) if path.isdir(data_dir / d)]
sample_dirs_empty = [data_dir_empty / d for d in listdir(data_dir_empty) if path.isdir(data_dir_empty / d)]


@pytest.fixture
def setup() -> Dict[str, Any]:
    return {
        'reader': SnippyVariantsReader(sample_dirs)
    }


def test_read_vcf(setup):
    vcf_file = data_dir / 'SampleA' / 'snps.vcf.gz'

    df = setup['reader'].read_vcf(vcf_file, 'SampleA')

    assert 46 == len(df), 'Data fram has incorrect length'

    assert {'snps.vcf.gz'} == set(df['FILE'].tolist()), 'Incorrect filename'
    assert {'SampleA'} == set(df['SAMPLE'].tolist()), 'Incorrect sample name'

    v = df[df['POS'] == 1048]
    assert 'C' == v['REF'].values[0], 'Incorrect reference'
    assert 'G' == v['ALT'].values[0], 'Incorrect alt'

    v = df[df['POS'] == 1253]
    assert 'T' == v['REF'].values[0], 'Incorrect reference'
    assert 'TAA' == v['ALT'].values[0], 'Incorrect alt'

    v = df[df['POS'] == 3656]
    assert 'CATT' == v['REF'].values[0], 'Incorrect reference'
    assert 'C' == v['ALT'].values[0], 'Incorrect alt'


def test_get_variants_table(setup):
    df = setup['reader'].get_features_table()

    assert 129 == len(df), 'Data has incorrect length'
    assert {'SampleA', 'SampleB', 'SampleC'} == set(df['SAMPLE'].tolist()), 'Incorrect sample names'


def test_read_genomic_masks_from_file(setup):
    sequence_file = data_dir / 'SampleA' / 'snps.aligned.fa'

    genomic_mask = setup['reader'].read_genomic_masks_from_file(sequence_file)

    assert 437 == len(genomic_mask)


def test_get_genomic_masks(setup):
    genomic_masks = setup['reader'].get_genomic_masked_regions()

    assert {'SampleA', 'SampleB', 'SampleC'} == set(genomic_masks.keys()), 'Incorrect samples'
    assert 437 == len(genomic_masks['SampleA'])
    assert 276 == len(genomic_masks['SampleB'])
    assert 329 == len(genomic_masks['SampleC'])


def test_get_samples_list(setup):
    reader = setup['reader']

    assert {'SampleA', 'SampleB', 'SampleC'} == set(reader.samples_list())


def test_get_samples_list_two_files():
    sample_dirs = [data_dir / 'SampleA', data_dir / 'SampleB']
    reader = SnippyVariantsReader(sample_dirs)

    assert {'SampleA', 'SampleB'} == set(reader.samples_list())


def test_read_empty_vcf():
    reader = SnippyVariantsReader(sample_dirs_empty)
    df = reader.get_features_table()

    assert 0 == len(df), 'Data has incorrect length'
