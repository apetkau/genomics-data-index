import warnings
from typing import Union, Set

import pytest
from pybedtools import BedTool

warnings.filterwarnings("ignore", category=DeprecationWarning)

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

from genomics_data_index.test.integration import data_dir, masked_positions_snippy_bed
from genomics_data_index.storage.service.CoreAlignmentService import CoreAlignmentService


def remove_column(alignment: MultipleSeqAlignment, position: int) -> MultipleSeqAlignment:
    if position <= 0 or position > alignment.get_alignment_length():
        raise Exception(f'position [{position}] is in an invalid range')

    return alignment[:, :(position - 1)] + alignment[:, position:]


def replace_column_with_reference(alignment: MultipleSeqAlignment, position: Union[int, range],
                                  skip_missing=True) -> MultipleSeqAlignment:
    if isinstance(position, range):
        for pos in position:
            alignment = replace_column_with_reference_internal(alignment, position=pos, skip_missing=skip_missing)
        return alignment
    else:
        return replace_column_with_reference_internal(alignment, position=position, skip_missing=skip_missing)


def replace_column_with_reference_internal(alignment: MultipleSeqAlignment, position: int,
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


def replace_gap_with_n_and_upper(alignment: MultipleSeqAlignment,
                                 position_set: Set[int] = None) -> MultipleSeqAlignment:
    for record in alignment:
        seq = record.seq.tomutable()
        for index, character in enumerate(seq):
            if position_set is None or (index + 1) in position_set:
                if character.upper() == '-' or character.upper() == 'X':
                    seq[index] = 'N'
                else:
                    seq[index] = seq[index].upper()

        record.seq = seq.toseq()
    return alignment


def replace_n_with_gap_and_upper(alignment: MultipleSeqAlignment,
                                 position_set: Set[int] = None) -> MultipleSeqAlignment:
    for record in alignment:
        seq = record.seq.tomutable()
        for index, character in enumerate(seq):
            if position_set is None or (index + 1) in position_set:
                if character.upper() == 'N':
                    seq[index] = '-'
                else:
                    seq[index] = seq[index].upper()

        record.seq = seq.toseq()
    return alignment


def alignment_upper(alignment: MultipleSeqAlignment):
    for record in alignment:
        record.seq = record.seq.upper()
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
def expected_alignment_full_unmodified() -> MultipleSeqAlignment:
    """
    Loads the expected alignment for testing (full alignment).
    This expected alignment was copied from the output of snippy.
    :return: The exepcted alignment as a MultipleSeqAlignment object.
    """
    with open(data_dir / 'snippy.full.aln', 'r') as f:
        alignments = list(AlignIO.parse(f, 'fasta'))
        alignment = alignments[0]

        replace_ref_name(alignment)
    return alignment


@pytest.fixture
def expected_alignment_full_only_snps(expected_alignment_full_unmodified) -> MultipleSeqAlignment:
    """
    Loads the full snippy alignment, and then removes all but SNPs.
    :return: The exepcted alignment as a MultipleSeqAlignment object.
    """
    alignment = expected_alignment_full_unmodified

    # Replace deletions/complex
    # Generated using
    # zgrep -hv '^#' Sample*/snps.vcf.gz | zgrep 'del\|complex' | cut -f 2,4 | sort -nu |
    #  sed -e 's/^\([0-9][0-9]*\)\t\(\w*\)/alignment =
    #   replace_column_with_reference(alignment, range(\1 + 1, \1 + 1 + len("\2") - 1), skip_missing=False)/'
    alignment = replace_column_with_reference(alignment, range(302 + 1, 302 + 1 + len("CT") - 1),
                                              skip_missing=False)
    alignment = replace_column_with_reference(alignment, range(347 + 1, 347 + 1 + len("TGAAG") - 1),
                                              skip_missing=False)
    alignment = replace_column_with_reference(alignment, range(349 + 1, 349 + 1 + len("AAGT") - 1),
                                              skip_missing=False)

    # This is a complex/OTHER event, which won't be represented in my alignment
    # Needed to expand boundaries on the left/right by 1 since it's not a simple deletion
    alignment = replace_column_with_reference(alignment, range(461, 461 + len("AAAT")),
                                              skip_missing=False)

    alignment = replace_column_with_reference(alignment, range(866 + 1, 866 + 1 + len("GCCAGATCC") - 1),
                                              skip_missing=False)
    alignment = replace_column_with_reference(alignment, range(883 + 1, 883 + 1 + len("CACATG") - 1),
                                              skip_missing=False)
    alignment = replace_column_with_reference(alignment, range(1070 + 1, 1070 + 1 + len("TG") - 1),
                                              skip_missing=False)
    alignment = replace_column_with_reference(alignment, range(1135 + 1, 1135 + 1 + len("CCT") - 1),
                                              skip_missing=False)
    alignment = replace_column_with_reference(alignment, range(1167 + 1, 1167 + 1 + len("GA") - 1),
                                              skip_missing=False)
    alignment = replace_column_with_reference(alignment, range(1209 + 1, 1209 + 1 + len("AG") - 1),
                                              skip_missing=False)
    alignment = replace_column_with_reference(alignment, range(1270 + 1, 1270 + 1 + len("CTT") - 1),
                                              skip_missing=False)

    # This is a complex/OTHER event, so won't be represented in my alignment. Need to treat left/right
    # boundaries differently
    alignment = replace_column_with_reference(alignment, range(1325, 1325 + len("TGA")),
                                              skip_missing=False)

    alignment = replace_column_with_reference(alignment,
                                              range(1483 + 1, 1483 + 1 + len("AAAGAGGGGCTGCTGGAGCCG") - 1),
                                              skip_missing=False)
    alignment = replace_column_with_reference(alignment, range(1708 + 1, 1708 + 1 + len("ATGCTGTTCAATAC") - 1),
                                              skip_missing=False)
    alignment = replace_column_with_reference(alignment, range(1815 + 1, 1815 + 1 + len("CAA") - 1),
                                              skip_missing=False)

    # This is a complex/OTHER event, so won't be represented in my alignment. Need to treat left/right
    # boundaries differently
    alignment = replace_column_with_reference(alignment, range(1984, 1984 + len("GTGATTG")),
                                              skip_missing=False)

    alignment = replace_column_with_reference(alignment, range(2419 + 1, 2419 + 1 + len("GA") - 1),
                                              skip_missing=False)
    alignment = replace_column_with_reference(alignment, range(2480 + 1, 2480 + 1 + len("TG") - 1),
                                              skip_missing=False)
    alignment = replace_column_with_reference(alignment, range(3008 + 1, 3008 + 1 + len("GA") - 1),
                                              skip_missing=False)
    alignment = replace_column_with_reference(alignment, range(3276 + 1, 3276 + 1 + len("GC") - 1),
                                              skip_missing=False)
    alignment = replace_column_with_reference(alignment, range(3576 + 1, 3576 + 1 + len("TTTC") - 1),
                                              skip_missing=False)
    alignment = replace_column_with_reference(alignment, range(3656 + 1, 3656 + 1 + len("CATT") - 1),
                                              skip_missing=False)
    alignment = replace_column_with_reference(alignment, range(3897 + 1, 3897 + 1 + len("GCGCA") - 1),
                                              skip_missing=False)
    alignment = replace_column_with_reference(alignment, range(4437 + 1, 4437 + 1 + len("AG") - 1),
                                              skip_missing=False)

    alignment = replace_gap_with_n_and_upper(alignment)

    alignment.sort()

    return alignment


@pytest.fixture
def expected_alignment_full_snps_mnp_deletions(expected_alignment_full_unmodified) -> MultipleSeqAlignment:
    """
    Loads the full snippy alignment, and then keeps SNPs, MNPs (which aren't in test data) and DELETIONs
    :return: The exepcted alignment as a MultipleSeqAlignment object.
    """
    alignment = alignment_upper(expected_alignment_full_unmodified)

    masked_positions = BedTool(masked_positions_snippy_bed)
    positions_to_fill_with_n = set()
    for interval in masked_positions:
        positions_to_fill_with_n = positions_to_fill_with_n | set(range(interval.start + 1, interval.end + 1))

    alignment = replace_gap_with_n_and_upper(alignment, position_set=positions_to_fill_with_n)

    # This is a complex/OTHER event, which won't be represented in my alignment
    # Needed to expand boundaries on the left/right by 1 since it's not a simple deletion
    alignment = replace_column_with_reference(alignment, range(461, 461 + len("AAAT")),
                                              skip_missing=False)

    # This is a complex/OTHER event, so won't be represented in my alignment. Need to treat left/right
    # boundaries differently
    alignment = replace_column_with_reference(alignment, range(1325, 1325 + len("TGA")),
                                              skip_missing=False)

    # This is a complex/OTHER event, so won't be represented in my alignment. Need to treat left/right
    # boundaries differently
    alignment = replace_column_with_reference(alignment, range(1984, 1984 + len("GTGATTG")),
                                              skip_missing=False)

    # These are recorded as deletions in the VCF file, but snippy is setting some characters to N instead of -
    # Since I set them as - in my results I need to change them in the expected alignment
    alignment = replace_n_with_gap_and_upper(alignment, {1210, 2420, 3657, 3658, 3659, 4438,
                                                         885, 1272, 2481,
                                                         1136, 1816, 3009})

    alignment.sort()

    return alignment


def get_unequal_positions_str(expected, actual) -> str:
    return_val = ''
    for idx, e in enumerate(expected.seq):
        a = actual[idx]
        if e != a:
            return_val = return_val + f'{idx + 1}: expected:{expected.id}={e} != actual:{actual.id}={a}\n'

    return return_val


def compare_alignments(expected: MultipleSeqAlignment, actual: MultipleSeqAlignment) -> None:
    assert len(expected) == len(actual), 'Alignment has incorrect number of samples'
    assert expected.get_alignment_length() == actual.get_alignment_length()

    number_of_sequences = len(expected)
    for i in range(0, number_of_sequences):
        assert expected[i].id == actual[i].id, f'Alignment ids are not equal [{expected[i].id}] != [{actual[i].id}]'
        assert expected[i].seq == actual[i].seq, f'Alignment sequences are not equal for ' \
                                                 f'[expected.id={expected[i].id}, actual.id={actual[i].id}] ' \
                                                 f'\n{get_unequal_positions_str(expected[i], actual[i])}'


def test_snippy_core_align(core_alignment_service: CoreAlignmentService, expected_alignment_core):
    actual_alignment = core_alignment_service.construct_alignment(reference_name='genome',
                                                                  samples=['SampleA', 'SampleB', 'SampleC'],
                                                                  include_variants=['SNP'])
    compare_alignments(expected_alignment_core, actual_alignment)


def test_snippy_core_align_with_other_include(core_alignment_service: CoreAlignmentService):
    with pytest.raises(Exception) as execinfo:
        core_alignment_service.construct_alignment(reference_name='genome',
                                                   samples=['SampleA', 'SampleB', 'SampleC'],
                                                   include_variants=['MNP'])
    assert 'Currently align_type=core only works with include_variants=["SNP"]' in str(execinfo.value)


def test_snippy_full_align_with_invalid_include_variants(core_alignment_service: CoreAlignmentService):
    with pytest.raises(Exception) as execinfo:
        core_alignment_service.construct_alignment(reference_name='genome',
                                                   samples=['SampleA', 'SampleB', 'SampleC'],
                                                   align_type='full',
                                                   include_variants=['invalid'])
    assert f"Only {['SNP', 'MNP', 'DELETION', 'DELETION_OTHER']} are supported" in str(execinfo.value)


def test_snippy_full_align(core_alignment_service, expected_alignment_full_only_snps):
    actual_alignment = core_alignment_service.construct_alignment(reference_name='genome',
                                                                  samples=['SampleA', 'SampleB', 'SampleC'],
                                                                  align_type='full',
                                                                  include_variants=['SNP'])
    compare_alignments(expected_alignment_full_only_snps, actual_alignment)


def test_snippy_full_align_snp_mnp_deletion(core_alignment_service, expected_alignment_full_snps_mnp_deletions):
    actual_alignment = core_alignment_service.construct_alignment(reference_name='genome',
                                                                  samples=['SampleA', 'SampleB', 'SampleC'],
                                                                  align_type='full')
    compare_alignments(expected_alignment_full_snps_mnp_deletions, actual_alignment)


def test_regular_vcf_full_align(core_alignment_service_non_snippy_vcfs, expected_alignment_full_only_snps):
    actual_alignment = core_alignment_service_non_snippy_vcfs.construct_alignment(reference_name='genome',
                                                                                  samples=['SampleA', 'SampleB',
                                                                                           'SampleC'],
                                                                                  align_type='full',
                                                                                  include_variants=['SNP'])
    compare_alignments(expected_alignment_full_only_snps, actual_alignment)


def test_regular_vcf_core_align(core_alignment_service_non_snippy_vcfs, expected_alignment_core):
    actual_alignment = core_alignment_service_non_snippy_vcfs.construct_alignment(reference_name='genome',
                                                                                  samples=['SampleA', 'SampleB',
                                                                                           'SampleC'],
                                                                                  include_variants=['SNP'])
    compare_alignments(expected_alignment_core, actual_alignment)
