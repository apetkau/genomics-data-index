from __future__ import annotations

import Bio.Seq
from bitarray import bitarray


class CoreBitMask:

    def __init__(self, sequence: Bio.Seq.Seq = None, existing_bitmask: bitarray = None):
        if existing_bitmask is not None and sequence is not None:
            raise Exception(f'Cannot set both existing_bitmask={existing_bitmask} and sequence={sequence}')

        if existing_bitmask:
            self._core_bitmask = existing_bitmask
        elif sequence:
            self._core_bitmask = bitarray(len(sequence))
            self._core_bitmask.setall(True)
            self._add_sequence(sequence)
        else:
            raise Exception('If no existing_bitmask set then sequence must be defined')

    def append_bitmask(self, bitmask: CoreBitMask) -> CoreBitMask:
        if bitmask is None:
            return self
        else:
            combined_array = self._core_bitmask & bitmask._core_bitmask
            return CoreBitMask(existing_bitmask=combined_array)

    def _add_sequence(self, sequence: Bio.Seq.Seq) -> None:
        for idx, char in enumerate(sequence):
            if char.upper() == 'N' or char == '-':
                self._core_bitmask[idx] = False

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
