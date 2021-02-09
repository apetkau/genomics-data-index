import warnings

import pytest

warnings.filterwarnings("ignore", category=DeprecationWarning)

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

from storage.test.integration.variant import data_dir


def remove_column(alignment: MultipleSeqAlignment, position: int) -> MultipleSeqAlignment:
    if position <= 0 or position > alignment.get_alignment_length():
        raise Exception(f'position [{position}] is in an invalid range')

    return alignment[:, :(position - 1)] + alignment[:, position:]


def replace_column_with_reference(alignment: MultipleSeqAlignment, position: int,
                                  skip_missing=True) -> MultipleSeqAlignment:
    column = alignment[:, (position - 1):position]
    ref = None
    ref_name = 'reference'
    for s in column:
        if s.id == ref_name:
            ref = s

    if ref is None:
        raise Exception('Could not find reference')

    for s in column:
        if s.id != ref_name:
            ss = s.seq.upper()
            if (not skip_missing) or (ss != 'N' and ss != '-'):
                s.seq = ref.seq

    return alignment[:, :(position - 1)] + column + alignment[:, position:]


def replace_ref_name(alignment):
    # Change reference name
    for s in alignment:
        s.description = 'generated automatically'
        if s.id == 'Reference':
            s.id = 'reference'

@pytest.fixture
def expected_alignment_core() -> MultipleSeqAlignment:
    '''
    Loads the expected alignment for testing (core alignment).
    This expected alignment was copied from the output of snippy. However, snippy
    includes some complex events in the alignment which I do not include.
    I modify this expected alignment (from snippy) by subtracting 2 SNVs/SNPs introduced
    by complex events. That way I can compare the alignments for equality.
    :return: The exepcted alignment as a MultipleSeqAlignment object.
    '''
    with open(data_dir / 'snippy.core.aln', 'r') as f:
        alignments = list(AlignIO.parse(f, 'fasta'))
        alignment = alignments[0]

        replace_ref_name(alignment)

        # Remove SNVs/SNPs introduced by complex events
        alignment = remove_column(alignment, 15)  # Position 1325
        alignment = remove_column(alignment, 23)  # Position 1984

        alignment.sort()

        return alignment


@pytest.fixture
def expected_alignment_full() -> MultipleSeqAlignment:
    '''
    Loads the expected alignment for testing (full alignment).
    This expected alignment was copied from the output of snippy. However, snippy
    includes some complex events in the alignment which I do not include.
    I modify this expected alignment (from snippy) by replacing 2 SNVs/SNPs introduced
    by complex events with the reference sequence. That way I can compare the alignments for equality.
    Additionally, snippy uses 'N' and '-' to mask regions, but I only use one character (default 'N')
    so I replace '-' with 'N'.
    :return: The exepcted alignment as a MultipleSeqAlignment object.
    '''
    with open(data_dir / 'snippy.full.aln', 'r') as f:
        alignments = list(AlignIO.parse(f, 'fasta'))
        alignment = alignments[0]

        replace_ref_name(alignment)

        # Replace SNVs/SNPs introduced by complex events
        alignment = replace_column_with_reference(alignment, 1325)
        alignment = replace_column_with_reference(alignment, 1984)

        alignment.sort()

        return alignment


def compare_alignments(a: MultipleSeqAlignment, b: MultipleSeqAlignment) -> None:
    assert len(a) == len(b), 'Alignment has incorrect number of samples'
    assert a.get_alignment_length() == b.get_alignment_length()

    alignment_length = a.get_alignment_length()
    for i in range(0, alignment_length):
        assert a[0].id == b[0].id, 'Alignment ids are not equal'
        assert a[0].seq == b[0].seq, 'Alignment sequences are not equal'


def test_snippy_align(core_alignment_service, expected_alignment_core):
    actual_alignment = core_alignment_service.construct_alignment(reference_name='genome',
                                                                  samples=['SampleA', 'SampleB', 'SampleC'])
    compare_alignments(expected_alignment_core, actual_alignment)


# def test_snippy_full_align(core_alignment_service, expected_alignment_full):
#     actual_alignment = core_alignment_service.construct_alignment(reference_name='genome',
#                                                                   samples=['SampleA', 'SampleB', 'SampleC'],
#                                                                   align_type='full')
#     compare_alignments(expected_alignment_full, actual_alignment)


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
