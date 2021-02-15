import math

from Bio.Seq import Seq

from storage.variant.CoreBitMask import CoreBitMask


def test_create_from_sequence():
    sequence = Seq('ATCG-NN')
    mask = CoreBitMask.from_sequence(sequence=sequence)

    assert 7 == len(mask), 'Invalid length'
    assert 4 == mask.core_length(), 'Invalid core length'
    assert math.isclose((4 / 7), mask.core_proportion()), 'Invalid core proportion'
    assert 1 in mask
    assert 4 in mask
    assert 5 not in mask
    assert 7 not in mask


def test_create_from_sequence2():
    sequence = Seq('---------')
    mask = CoreBitMask.from_sequence(sequence=sequence)

    assert 9 == len(mask), 'Invalid length'
    assert 0 == mask.core_length(), 'Invalid core length'
    assert 0 == mask.core_proportion()
    assert 1 not in mask
    assert 9 not in mask


def test_append():
    sequence1 = Seq('ATCG-NN')
    mask1 = CoreBitMask.from_sequence(sequence=sequence1)
    sequence2 = Seq('NTCGATT')
    mask2 = CoreBitMask.from_sequence(sequence=sequence2)

    assert 7 == len(mask1), 'Invalid length'
    assert 4 == mask1.core_length(), 'Invalid core length'
    assert 7 == len(mask2), 'Invalid length'
    assert 6 == mask2.core_length(), 'Invalid core length'

    combined = mask1.append_bitmask(mask2)
    assert 7 == len(combined)
    assert 3 == combined.core_length()

    assert math.isclose((3 / 7), combined.core_proportion())

    assert 1 not in combined
    assert 2 in combined
    assert 4 in combined
    assert 5 not in combined
    assert 7 not in combined


def test_iter_missing():
    sequence = Seq('ATNG-NN')
    mask = CoreBitMask.from_sequence(sequence)

    assert [3,5,6,7] == list(mask.iter_missing_positions())
