import pytest
import warnings

warnings.filterwarnings("ignore", category=DeprecationWarning)

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

from storage.test.integration.variant import data_dir


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
        alignment = remove_column(alignment, 15)  # Position 1325
        alignment = remove_column(alignment, 23)  # Position 1984

        alignment.sort()

        return alignment


def compare_alignments(a: MultipleSeqAlignment, b: MultipleSeqAlignment) -> None:
    assert len(a) == len(b), 'Alignment has incorrect number of samples'
    assert a.get_alignment_length() == b.get_alignment_length()

    alignment_length = a.get_alignment_length()
    for i in range(0, alignment_length):
        assert a[0].id == b[0].id, 'Alignment ids are not equal'
        assert a[0].seq == b[0].seq, 'Alignment sequences are not equal'


def test_snippy_align(core_alignment_service, expected_alignment):
    actual_alignment = core_alignment_service.construct_alignment(reference_name='genome',
                                                                  samples=['SampleA', 'SampleB', 'SampleC'])
    compare_alignments(expected_alignment, actual_alignment)


def test_get_variants(core_alignment_service):
    variants = core_alignment_service._get_variants(sequence_name='reference')

    assert 60 == len(variants.keys()), 'Incorrect number of variants returned (counting only SNV/SNPs)'

    assert 'reference:4265:G:C' == variants[4265]['SampleA'].id, 'Incorrect variant returned'
    assert 'SampleB' not in variants[4265], 'Incorrect variant returned'
    assert 'SampleC' not in variants[4265], 'Incorrect variant returned'

    assert 'reference:839:C:G' == variants[839]['SampleB'].id, 'Incorrect variant returned'
    assert 'reference:839:C:G' == variants[839]['SampleC'].id, 'Incorrect variant returned'
    assert 'SampleA' not in variants[839], 'Incorrect variant returned'


def test_sample_sequence(core_alignment_service):
    sample_sequences = core_alignment_service._sample_sequence(reference_name='genome',
                                                               samples=['SampleA', 'SampleB', 'SampleC'])

    assert {'reference'} == set(sample_sequences.keys()), 'Incorrect number of sequences'

    sample_names = {ss.sample.name for ss in sample_sequences['reference']}
    assert {'SampleA', 'SampleB', 'SampleC'} == sample_names
