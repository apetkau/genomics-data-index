import warnings

import pytest

warnings.filterwarnings("ignore", category=DeprecationWarning)

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

from genomics_data_index.test.integration import data_dir
from genomics_data_index.storage.service.CoreAlignmentService import CoreAlignmentService


def remove_column(alignment: MultipleSeqAlignment, position: int) -> MultipleSeqAlignment:
    if position <= 0 or position > alignment.get_alignment_length():
        raise Exception(f'position [{position}] is in an invalid range')

    return alignment[:, :(position - 1)] + alignment[:, position:]


def replace_column_with_reference(alignment: MultipleSeqAlignment, position: int,
                                  skip_missing=True) -> MultipleSeqAlignment:
    column = alignment[:, (position - 1):position]
    ref = None
    ref_name = 'genome'
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
        if s.id == 'Reference':
            s.id = 'genome'
            s.description = '[reference genome]'
        else:
            s.description = 'generated automatically'


def replace_gap_with_n_and_upper(alignment: MultipleSeqAlignment) -> MultipleSeqAlignment:
    for record in alignment:
        seq = record.seq.tomutable()
        for index, character in enumerate(seq):
            if character.upper() == '-' or character.upper() == 'X':
                seq[index] = 'N'
            else:
                seq[index] = seq[index].upper()

        record.seq = seq.toseq()
    return alignment


@pytest.fixture
def expected_alignment_core() -> MultipleSeqAlignment:
    '''
    Loads the expected alignment for testing (core alignment).
    This expected alignment was copied from the output of snippy. However, snippy
    includes some complex events in the alignment which I do not include.
    I modify this expected alignment (from snippy) by subtracting 2 SNVs/SNPs introduced
    by complex events. That way I can compare the alignments for equality.
    :return: The expected alignment as a MultipleSeqAlignment object.
    '''
    with open(data_dir / 'snippy.core.aln', 'r') as f:
        alignments = list(AlignIO.parse(f, 'fasta'))
        alignment = alignments[0]

        replace_ref_name(alignment)

        # Remove SNVs/SNPs introduced by complex events
        alignment = remove_column(alignment, 15)  # Position 1325
        alignment = remove_column(alignment, 23 - 1)  # Position 1984 (-1 to account for previous removed column)

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

        # Replace Ns in locations I did not expect (not low coverage)
        alignment = replace_column_with_reference(alignment, 463, skip_missing=False)
        alignment = replace_column_with_reference(alignment, 464, skip_missing=False)

        # Replace SNVs/SNPs introduced by complex events
        alignment = replace_column_with_reference(alignment, 1325)
        alignment = replace_column_with_reference(alignment, 1984)
        alignment = replace_gap_with_n_and_upper(alignment)

        alignment.sort()

        return alignment


def compare_alignments(a: MultipleSeqAlignment, b: MultipleSeqAlignment) -> None:
    assert len(a) == len(b), 'Alignment has incorrect number of samples'
    assert a.get_alignment_length() == b.get_alignment_length()

    number_of_sequences = len(a)
    for i in range(0, number_of_sequences):
        assert a[i].id == b[i].id, f'Alignment ids are not equal [{a[i].id}] != [{b[i].id}]'
        assert a[i].seq == b[i].seq, f'Alignment sequences are not equal for [a.id={a[i].id}, b.id={b[i].id}]'


def test_snippy_core_align(core_alignment_service: CoreAlignmentService, expected_alignment_core):
    actual_alignment = core_alignment_service.construct_alignment(reference_name='genome',
                                                                  samples=['SampleA', 'SampleB', 'SampleC'])
    compare_alignments(expected_alignment_core, actual_alignment)


def test_snippy_full_align(core_alignment_service, expected_alignment_full):
    actual_alignment = core_alignment_service.construct_alignment(reference_name='genome',
                                                                  samples=['SampleA', 'SampleB', 'SampleC'],
                                                                  align_type='full')
    compare_alignments(expected_alignment_full, actual_alignment)


def test_regular_vcf_full_align(core_alignment_service_non_snippy_vcfs, expected_alignment_full):
    actual_alignment = core_alignment_service_non_snippy_vcfs.construct_alignment(reference_name='genome',
                                                                                  samples=['SampleA', 'SampleB',
                                                                                           'SampleC'],
                                                                                  align_type='full')
    compare_alignments(expected_alignment_full, actual_alignment)


def test_regular_vcf_core_align(core_alignment_service_non_snippy_vcfs, expected_alignment_core):
    actual_alignment = core_alignment_service_non_snippy_vcfs.construct_alignment(reference_name='genome',
                                                                                  samples=['SampleA', 'SampleB',
                                                                                           'SampleC'])
    compare_alignments(expected_alignment_core, actual_alignment)
