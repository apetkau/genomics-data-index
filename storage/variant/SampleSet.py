from __future__ import annotations

from typing import Iterable, Generator, Union, Set

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

    def intersection(self, other: Union[Set[int]]) -> SampleSet:
        if other is None:
            raise Exception('Cannot intersect [other = None]')
        elif isinstance(other, SampleSet):
            return self._bitmap.intersection(other._bitmap)
        elif isinstance(other, set):
            return self._bitmap.intersection(BitMap(other))
        else:
            raise Exception(f'Cannot intersect other of type [{type(other)}]')

    @staticmethod
    def from_bytes(data: bytes) -> SampleSet:
        bitmap = BitMap.deserialize(data)
        return SampleSet(existing_bitmap=bitmap)

    def get_bytes(self) -> bytes:
        return self._bitmap.serialize()

    def __iter__(self) -> Generator[int, None, None]:
        yield from self._bitmap

    def __contains__(self, value: int) -> bool:
        return value in self._bitmap

    def __len__(self) -> int:
        return len(self._bitmap)

    def __repr__(self):
        return (f'<SampleSet(size={len(self._bitmap)}>')
