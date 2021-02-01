import tempfile
import warnings
from pathlib import Path

warnings.filterwarnings("ignore", category=DeprecationWarning)

import pytest
from Bio import AlignIO

from storage.test.integration.variant import data_dir, sample_dirs, reference_file
from storage.variant.service import DatabaseConnection
from storage.variant.service.CoreAlignmentService import CoreAlignmentService
from storage.variant.service.ReferenceService import ReferenceService
from storage.variant.service.VariationService import VariationService
from storage.variant.VariantsReader import SnippyVariantsReader


@pytest.fixture
def database() -> DatabaseConnection:
    return DatabaseConnection('sqlite:///:memory:')


@pytest.fixture
def reference_service(database) -> ReferenceService:
    seq_repo_root = Path(tempfile.mkdtemp(prefix='index-test'))
    reference_service = ReferenceService(database, seq_repo_root)

    reference_service.add_reference_genome(reference_file)

    return reference_service


@pytest.fixture
def variation_service(database, reference_service) -> VariationService:
    reference_service.create_reference_genome(reference_file)
    reference = reference_service.find_reference_genome('genome')

    ref_contigs = {s.sequence_name: s for s in reference.sequences}

    variants_reader = SnippyVariantsReader(sample_dirs)

    var_df = variants_reader.get_variants_table()
    core_masks = variants_reader.get_core_masks()

    var_service = VariationService(database)
    var_service.insert_variants(var_df=var_df,
                                ref_contigs=ref_contigs,
                                core_masks=core_masks)

    return var_service


@pytest.fixture
def core_alignment_service(database, reference_service) -> CoreAlignmentService:
    return CoreAlignmentService(database, reference_service)


def test_snippy_align(core_alignment_service, variation_service):
    print(variation_service)

    actual_alignment = core_alignment_service.construct_alignment(reference_name='genome',
                                                                  samples=['SampleA', 'SampleB', 'SampleC'])

    with open(data_dir / 'snippy.core.aln', 'r') as f:
        alignments = list(AlignIO.parse(f, 'fasta'))
        expected_alignment = alignments[0]

    assert len(expected_alignment) == len(actual_alignment), 'Alignment has incorrect length'


def test_get_variants(core_alignment_service, variation_service):
    variants = core_alignment_service._get_variants(sequence_name='reference')

    assert 60 == len(variants.keys()), 'Incorrect number of variants returned (counting only SNV/SNPs)'

    assert 'reference:4265:G:C' == variants[4265]['SampleA'].id, 'Incorrect variant returned'
    assert 'SampleB' not in variants[4265], 'Incorrect variant returned'
    assert 'SampleC' not in variants[4265], 'Incorrect variant returned'

    assert 'reference:839:C:G' == variants[839]['SampleB'].id, 'Incorrect variant returned'
    assert 'reference:839:C:G' == variants[839]['SampleC'].id, 'Incorrect variant returned'
    assert 'SampleA' not in variants[839], 'Incorrect variant returned'


def test_sample_sequence(core_alignment_service, variation_service):
    sample_sequences = core_alignment_service._sample_sequence(reference_name='genome',
                                                               samples=['SampleA', 'SampleB', 'SampleC'])

    assert {'reference'} == set(sample_sequences.keys()), 'Incorrect number of sequences'

    sample_names = {ss.sample.name for ss in sample_sequences['reference']}
    assert {'SampleA', 'SampleB', 'SampleC'} == sample_names
