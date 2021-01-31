from typing import Dict, Any
import pytest

from os import path, listdir
from pathlib import Path

from storage.variant.VariantsReader import SnippyVariantsReader

data_dir = Path(path.dirname(__file__), '..', 'data', 'snippy')
sample_dirs = [data_dir / d for d in listdir(data_dir) if path.isdir(data_dir / d)]


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
    df = setup['reader'].get_variants_table()

    assert 129 == len(df), 'Data has incorrect length'
    assert {'SampleA', 'SampleB', 'SampleC'} == set(df['SAMPLE'].tolist()), 'Incorrect sample names'


def test_read_core_masks(setup):
    sequence_file = data_dir / 'SampleA' / 'snps.aligned.fa'

    core_masks = setup['reader'].read_core_masks(sequence_file)

    total_length = len(core_masks['reference'])
    assert 5180 == total_length, f'File has incorrect total length [{total_length}]'

    missing_length = total_length - core_masks['reference'].core_length()
    assert 437 == missing_length, f'File has incorrect missing length [{missing_length}]'


def test_get_core_masks(setup):
    core_masks = setup['reader'].get_core_masks()

    assert {'SampleA', 'SampleB', 'SampleC'} == set(core_masks.keys()), 'Incorrect samples'
    assert 4743 == core_masks['SampleA']['reference'].core_length(), 'Incorrect core length for SampleA'
    assert 4904 == core_masks['SampleB']['reference'].core_length(), 'Incorrect core length for SampleB'
    assert 4851 == core_masks['SampleC']['reference'].core_length(), 'Incorrect core length for SampleC'