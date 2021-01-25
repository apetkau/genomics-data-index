from Bio.Seq import Seq
import math

from storage.variant.CoreBitMask import CoreBitMask


def test_core_length():
    sequence = Seq('ATCG-NN')
    mask = CoreBitMask(sequence=sequence)

    assert 4 == mask.core_length(), 'Invalid core length'


def test_core_proportion():
    sequence = Seq('ATCG-NN')
    mask = CoreBitMask(sequence=sequence)

    assert math.isclose((4/7),mask.core_proportion()), 'Invalid core proportion'