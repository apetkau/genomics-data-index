from __future__ import annotations

from typing import Iterable, Generator, Union, Set, List

from pyroaring import BitMap


class SampleSet:

    def __init__(self, sample_ids: Iterable[int] = None, existing_bitmap: BitMap = None):
        if sample_ids is None and existing_bitmap is None:
            raise Exception('Both sample_ids and existing_bitmap are unset. One of them must be set.')
        elif sample_ids is not None and existing_bitmap is not None:
            raise Exception('Both sample_ids and existing_bitmap are set. Only one of them can be set')
        elif sample_ids is not None:
            self._bitmap = BitMap(sample_ids)
        else:
            self._bitmap = existing_bitmap

    def intersection(self, other: Union[Set[int], SampleSet]) -> SampleSet:
        if other is None:
            raise Exception('Cannot intersect other[None]')
        elif isinstance(other, AllSampleSet):
            return self
        elif isinstance(other, SampleSet):
            return SampleSet(existing_bitmap=self._bitmap.intersection(other._bitmap))
        elif isinstance(other, set):
            return SampleSet(self._bitmap.intersection(BitMap(other)))
        else:
            raise Exception(f'Cannot intersect other of type [{type(other)}]')

    def union(self, other: Union[Set[int], SampleSet]) -> SampleSet:
        if other is None:
            raise Exception('Cannot union other=[None]')
        elif isinstance(other, AllSampleSet):
            return other
        elif isinstance(other, SampleSet):
            return SampleSet(existing_bitmap=self._bitmap.union(other._bitmap))
        elif isinstance(other, set):
            return SampleSet(self._bitmap.union(BitMap(other)))
        else:
            raise Exception(f'Cannot union other of type [{type(other)}]')

    def minus(self, other: Union[Set[int], SampleSet, List[SampleSet]]) -> SampleSet:
        if other is None:
            raise Exception('Cannot union other=[None]')
        elif isinstance(other, AllSampleSet):
            return SampleSet.create_empty()
        elif isinstance(other, SampleSet):
            return SampleSet(existing_bitmap=(self._bitmap - other._bitmap))
        elif isinstance(other, set):
            return SampleSet(existing_bitmap=(self._bitmap - BitMap(other)))
        elif isinstance(other, list):
            difference_bitmap = self._bitmap
            for s in other:
                if isinstance(s, SampleSet):
                    difference_bitmap = difference_bitmap - s._bitmap
                else:
                    raise Exception(f'Got invalid type=[{type(s)}] in input other=[{other}]. '
                                    f'Expected a list of SampleSet')
            return SampleSet(existing_bitmap=difference_bitmap)
        else:
            raise Exception(f'Cannot union other of type [{type(other)}]')

    def jaccard_index(self, other: Union[Set[int], SampleSet]) -> float:
        if other is None:
            raise Exception('Cannot calculate jaccard with other=[None]')
        elif isinstance(other, AllSampleSet):
            raise NotImplementedError('Jaccard with AllSampleSet is not implemented')
        elif isinstance(other, SampleSet):
            return self._bitmap.jaccard_index(other._bitmap)
        elif isinstance(other, set):
            return self._bitmap.jaccard_index(SampleSet(other)._bitmap)
        else:
            raise Exception(f'Cannot union other of type [{type(other)}]')

    def is_empty(self) -> bool:
        return len(self._bitmap) == 0

    @staticmethod
    def from_bytes(data: bytes) -> SampleSet:
        bitmap = BitMap.deserialize(data)
        return SampleSet(existing_bitmap=bitmap)

    @classmethod
    def create_empty(cls):
        return SampleSet(existing_bitmap=BitMap())

    @classmethod
    def create_all(cls):
        return AllSampleSet()

    def get_bytes(self) -> bytes:
        return self._bitmap.serialize()

    def __iter__(self) -> Generator[int, None, None]:
        yield from self._bitmap

    def __contains__(self, value: int) -> bool:
        return value in self._bitmap

    def __len__(self) -> int:
        return len(self._bitmap)

    def __repr__(self):
        return (f'<SampleSet(size={len(self._bitmap)})>')


class AllSampleSet(SampleSet):
    # Roaring bitmaps store 32-bit integers, so max size is all 32-bit ints
    MAX_VALUE = 2 ** 32

    def __init__(self):
        super().__init__(sample_ids=[])

    def intersection(self, other: Union[Set[int], SampleSet]) -> SampleSet:
        return other

    def is_empty(self) -> bool:
        return False

    def minus(self, other: Union[Set[int], SampleSet, List[SampleSet]]) -> SampleSet:
        raise Exception('Cannot subtract anything from AllSampleSet')

    def get_bytes(self) -> bytes:
        raise Exception('Cannot serialize AllSampleSet')

    def __iter__(self) -> Generator[int, None, None]:
        yield from range(self.MAX_VALUE)

    def __contains__(self, value: int) -> bool:
        return value < self.MAX_VALUE

    def __len__(self) -> int:
        return self.MAX_VALUE

    def __repr__(self):
        return (f'<AllSampleSet(size={len(self._bitmap)})>')
