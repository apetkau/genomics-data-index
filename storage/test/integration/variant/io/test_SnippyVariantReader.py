from os import path, listdir
from typing import Dict, Any
import tempfile
from pathlib import Path

import pytest

from storage.test.integration.variant import data_dir
from storage.test.integration.variant import data_dir_empty
from storage.variant.io.mutation.SnippyVariantsReader import SnippyVariantsReader
from storage.variant.io.mutation.VcfVariantsReader import VcfVariantsReader
from storage.variant.io.processor.SerialSampleFilesProcessor import SerialSampleFilesProcessor

sample_dirs = [data_dir / d for d in listdir(data_dir) if path.isdir(data_dir / d)]
sample_dirs_empty = [data_dir_empty / d for d in listdir(data_dir_empty) if path.isdir(data_dir_empty / d)]


@pytest.fixture
def variants_reader() -> SnippyVariantsReader:
    tmp_dir = Path(tempfile.mkdtemp())
    return SnippyVariantsReader.create(sample_dirs,
                                       sample_files_processor=SerialSampleFilesProcessor(tmp_dir))


def test_read_vcf(variants_reader):
    vcf_file = data_dir / 'SampleA' / 'snps.vcf.gz'

    df = variants_reader.read_vcf(vcf_file, 'SampleA')

    assert 46 == len(df), 'Data fram has incorrect length'

    assert {'snps.vcf.gz'} == set(df['FILE'].tolist()), 'Incorrect filename'
    assert {'SampleA'} == set(df['SAMPLE'].tolist()), 'Incorrect sample name'

    v = df[df['POS'] == 461]
    assert 'AAAT' == v['REF'].values[0], 'Incorrect reference'
    assert 'G' == v['ALT'].values[0], 'Incorrect alt'

    v = df[df['POS'] == 1048]
    assert 'C' == v['REF'].values[0], 'Incorrect reference'
    assert 'G' == v['ALT'].values[0], 'Incorrect alt'

    v = df[df['POS'] == 1253]
    assert 'T' == v['REF'].values[0], 'Incorrect reference'
    assert 'TAA' == v['ALT'].values[0], 'Incorrect alt'

    v = df[df['POS'] == 3656]
    assert 'CATT' == v['REF'].values[0], 'Incorrect reference'
    assert 'C' == v['ALT'].values[0], 'Incorrect alt'


def test_get_variants_table(variants_reader):
    df = variants_reader.get_features_table()

    assert 129 == len(df), 'Data has incorrect length'
    assert {'SampleA', 'SampleB', 'SampleC'} == set(df['SAMPLE'].tolist()), 'Incorrect sample names'


def test_get_genomic_masks(variants_reader):
    mask = variants_reader.get_genomic_masked_region('SampleA')
    assert 437 == len(mask)
    assert {'reference'} == mask.sequence_names()

    mask = variants_reader.get_genomic_masked_region('SampleB')
    assert 276 == len(mask)
    assert {'reference'} == mask.sequence_names()

    mask = variants_reader.get_genomic_masked_region('SampleC')
    assert 329 == len(mask)
    assert {'reference'} == mask.sequence_names()


def test_get_samples_list(variants_reader):
    assert {'SampleA', 'SampleB', 'SampleC'} == set(variants_reader.samples_list())


def test_get_samples_list_two_files():
    sample_dirs = [data_dir / 'SampleA', data_dir / 'SampleB']
    reader = SnippyVariantsReader.create(sample_dirs)

    assert {'SampleA', 'SampleB'} == set(reader.samples_list())


def test_read_empty_vcf():
    tmp_dir = Path(tempfile.mkdtemp())
    reader = SnippyVariantsReader.create(sample_dirs_empty,
                                         sample_files_processor=SerialSampleFilesProcessor(tmp_dir))
    df = reader.get_features_table()

    assert 0 == len(df), 'Data has incorrect length'
