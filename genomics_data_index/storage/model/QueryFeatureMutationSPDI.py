from __future__ import annotations

from typing import Union

from genomics_data_index.storage.model import NUCLEOTIDE_UNKNOWN
from genomics_data_index.storage.model.QueryFeature import QueryFeature


class QueryFeatureMutationSPDI(QueryFeature):

    def __init__(self, spdi: str):
        super().__init__()
        self._spdi = spdi

        seq, pos, ref, alt = self._spdi.split(self.SPLIT_CHAR)

        if seq == self.WILD:
            raise Exception(f'Cannot set seq to be wild ({self.WILD}): {spdi}')
        else:
            self._seq_name = seq

        if pos == self.WILD:
            raise Exception(f'Cannot set pos to be wild ({self.WILD}): {spdi}')
        else:
            self._pos = int(pos)

        if ref == self.WILD:
            raise Exception(f'Cannot set ref to be wild ({self.WILD}): {spdi}')
        elif ref.isdigit():
            self._ref = int(ref)
        else:
            self._ref = ref

        if alt is None:
            self._alt = self.WILD
        else:
            self._alt = alt

    @property
    def id(self) -> str:
        return self._spdi

    @property
    def scope(self) -> str:
        return self._seq_name

    @property
    def sequence(self) -> str:
        return self.scope

    @property
    def position(self) -> int:
        return self._pos

    @property
    def start(self) -> int:
        return self.position

    @property
    def stop(self) -> int:
        if isinstance(self._ref, int):
            return self.position + self._ref
        else:
            return self.position + len(self.ref)

    @property
    def start0(self) -> int:
        return self.start - 1

    @property
    def stop0(self) -> int:
        return self.stop - 1

    @property
    def ref(self) -> Union[str, int]:
        return self._ref

    @property
    def deletion(self) -> Union[str, int]:
        return self.ref

    @property
    def alt(self) -> str:
        return self._alt

    @property
    def insertion(self) -> str:
        return self.alt

    def is_unknown(self) -> bool:
        return False

    def to_unknown(self) -> QueryFeature:
        return QueryFeatureMutationSPDI(':'.join([
            self._seq_name,
            self._pos,
            self._ref,
            NUCLEOTIDE_UNKNOWN
        ]))
