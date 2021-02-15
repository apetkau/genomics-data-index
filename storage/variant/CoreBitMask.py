from __future__ import annotations

import Bio.Seq
from bitarray import bitarray


class CoreBitMask:

    def __init__(self, existing_bitmask: bitarray = None):
        if existing_bitmask:
            self._core_bitmask = existing_bitmask
        else:
            raise Exception('If no existing_bitmask set then sequence must be defined')

    def append_bitmask(self, bitmask: CoreBitMask) -> CoreBitMask:
        if bitmask is None:
            return self
        else:
            combined_array = self._core_bitmask & bitmask._core_bitmask
            return CoreBitMask(existing_bitmask=combined_array)

    @staticmethod
    def from_bytes(data: bytes, length: int) -> CoreBitMask:
        if length <= 0:
            raise Exception(f'Invalid length=[{length}]')
        else:
            barray = bitarray()
            barray.frombytes(data)

            # Since I'm decoding a bitarray from bytes, the bytes must always be a multiple of 8
            # But the sequence length can be a non-multiple of 8
            # So I have to remove (slice) the few additional elements from this array that got added
            # When decoding from bytes
            barray = barray[:length]
            return CoreBitMask(existing_bitmask=barray)

    @staticmethod
    def from_sequence(sequence: Bio.Seq.Seq):
        barray = bitarray(len(sequence))
        barray.setall(True)
        for idx, char in enumerate(sequence):
            if char.upper() == 'N' or char == '-':
                barray[idx] = False
        return CoreBitMask(existing_bitmask=barray)

    def get_bytes(self):
        return self._core_bitmask.tobytes()

    def core_length(self) -> int:
        return self._core_bitmask.count()

    def core_proportion(self) -> float:
        return self.core_length() / len(self)

    def __len__(self) -> int:
        return len(self._core_bitmask)

    def __getitem__(self, index: int) -> bool:
        return self._core_bitmask[index]
