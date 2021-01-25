from os import path
from pathlib import Path

from storage.variant.VariantsReader import SnippyVariantsReader

data_dir = Path(path.dirname(__file__), '..', 'data', 'snippy')


def test_read_vcf():
    reader = SnippyVariantsReader(data_dir)
    vcf_file = data_dir / 'SampleA' / 'snps.vcf.gz'

    df = reader.read_vcf(vcf_file)

    assert 46 == len(df), 'Data fram has incorrect length'


def test_read_core_masks():
    reader = SnippyVariantsReader(data_dir)
    sequence_file = data_dir / 'SampleA' / 'snps.aligned.fa'

    core_masks = reader.read_core_masks(sequence_file)

    total_length = len(core_masks['reference'])
    assert 5180 == total_length, f'File has incorrect total length [{total_length}]'

    missing_length = total_length - core_masks['reference'].core_length()
    assert 437 == missing_length, f'File has incorrect missing length [{missing_length}]'