import math

from Bio.Seq import Seq

from storage.variant.CoreBitMask import CoreBitMask


def test_create_from_sequence():
    sequence = Seq('ATCG-NN')
    mask = CoreBitMask.from_sequence(sequence=sequence)

    assert 7 == len(mask), 'Invalid length'
    assert 4 == mask.core_length(), 'Invalid core length'
    assert math.isclose((4 / 7), mask.core_proportion()), 'Invalid core proportion'


def test_create_from_sequence2():
    sequence = Seq('---------')
    mask = CoreBitMask.from_sequence(sequence=sequence)

    assert 9 == len(mask), 'Invalid length'
    assert 0 == mask.core_length(), 'Invalid core length'
    assert 0 == mask.core_proportion()


def test_create_from_bytes():
    mask = CoreBitMask.from_bytes(b'\xf0', 8)

    assert 8 == len(mask), 'Invalid length'
    assert 4 == mask.core_length(), 'Invalid core length'
    assert math.isclose((4 / 8), mask.core_proportion()), 'Invalid core proportion'


def test_create_from_bytes2():
    mask = CoreBitMask.from_bytes(b'\xf0', 7)

    assert 7 == len(mask), 'Invalid length'
    assert 4 == mask.core_length(), 'Invalid core length'
    assert math.isclose((4 / 7), mask.core_proportion()), 'Invalid core proportion'