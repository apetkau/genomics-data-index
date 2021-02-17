from __future__ import annotations
from typing import Generator

import Bio.Seq
from pyroaring import BitMap


class CoreBitMask:

    def __init__(self, existing_bitmask: BitMap, sequence_length: int):
        if existing_bitmask is None:
            raise Exception('Cannot create CoreBitMask with empty bitmask')
        elif not sequence_length or sequence_length <= 0:
            raise Exception(f'Invalid sequence_length=[{sequence_length}]')

        self._core_bitmask = existing_bitmask
        self._sequence_length = sequence_length

    def append_bitmask(self, bitmask: CoreBitMask) -> CoreBitMask:
        if bitmask is None:
            return self
        elif bitmask._sequence_length != self._sequence_length:
            raise Exception(f'Cannot append bitmask, sequence lengths'
                            f' [left={self._sequence_length}, right={bitmask._sequence_length}] are unequal')
        else:
            combined_array = self._core_bitmask | bitmask._core_bitmask
            return CoreBitMask(existing_bitmask=combined_array, sequence_length=self._sequence_length)

    def iter_missing_positions(self) -> Generator[int, None, None]:
        yield from self._core_bitmask

    @staticmethod
    def from_bytes(data: bytes, length: int) -> CoreBitMask:
        if length <= 0:
            raise Exception(f'Invalid length=[{length}]')
        else:
            bitmap = BitMap.deserialize(data)
            return CoreBitMask(existing_bitmask=bitmap, sequence_length=length)

    @staticmethod
    def from_sequence(sequence: Bio.Seq.Seq):
        bitmap = BitMap()
        for idx, char in enumerate(sequence):
            if char.upper() == 'N' or char == '-':
                # Add 1 so that positions are counted from 1 instead of 0
                bitmap.add(idx+1)
        return CoreBitMask(existing_bitmask=bitmap, sequence_length=len(sequence))

    @staticmethod
    def empty_mask(length: int):
        bitmap = BitMap()
        return CoreBitMask(existing_bitmask=bitmap, sequence_length=length)

    def get_bytes(self) -> bytes:
        return self._core_bitmask.serialize()

    def core_length(self) -> int:
        return self._sequence_length - len(self._core_bitmask)

    def is_empty(self):
        return len(self._core_bitmask) == 0

    def core_proportion(self) -> float:
        return self.core_length() / len(self)

    def __contains__(self, value: int) -> bool:
        return value not in self._core_bitmask

    def __len__(self) -> int:
        return self._sequence_length
