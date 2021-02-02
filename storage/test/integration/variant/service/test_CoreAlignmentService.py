import tempfile
import warnings
from pathlib import Path

warnings.filterwarnings("ignore", category=DeprecationWarning)

import pytest
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

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


@pytest.fixture
def expected_alignment() -> MultipleSeqAlignment:
    '''
    Loads the expected alignment for testing.
    This expected alignment was copied from the output of snippy. However, snippy
    includes some complex events in the alignment which I do not include.
    I modify this expected alignment (from snippy) by subtracting 2 SNVs/SNPs introduced
    by complex events. That way I can compare the alignments for equality.
    :return: The exepcted alignment as a MultipleSeqAlignment object.
    '''
    def remove_column(alignment: MultipleSeqAlignment, position: int) -> MultipleSeqAlignment:
        if position <= 0 or position > alignment.get_alignment_length():
            raise Exception(f'position [{position}] is in an invalid range')

        return alignment[:, :(position - 1)] + alignment[:, position:]

    with open(data_dir / 'snippy.core.aln', 'r') as f:
        alignments = list(AlignIO.parse(f, 'fasta'))
        alignment = alignments[0]

        # Change reference name
        for s in alignment:
            s.description = 'generated automatically'
            if s.id == 'Reference':
                s.id = 'reference'

        # Remove SNVs/SNPs introduced by complex events
        alignment = remove_column(alignment, 15) # Position 1325
        alignment = remove_column(alignment, 23) # Position 1984

        alignment.sort()

        return alignment


def compare_alignments(a: MultipleSeqAlignment, b: MultipleSeqAlignment) -> None:
    assert len(a) == len(b), 'Alignment has incorrect number of samples'
    assert a.get_alignment_length() == b.get_alignment_length()

    alignment_length = a.get_alignment_length()
    for i in range(0,alignment_length):
        assert a[0].id == b[0].id, 'Alignment ids are not equal'
        assert a[0].seq == b[0].seq, 'Alignment sequences are not equal'
        

def test_snippy_align(core_alignment_service, variation_service, expected_alignment):
    actual_alignment = core_alignment_service.construct_alignment(reference_name='genome',
                                                                  samples=['SampleA', 'SampleB', 'SampleC'])
    compare_alignments(expected_alignment, actual_alignment)


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
