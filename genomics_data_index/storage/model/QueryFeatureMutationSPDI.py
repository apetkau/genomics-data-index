from __future__ import annotations

from typing import Union, List, Optional

from genomics_data_index.storage.model import NUCLEOTIDE_UNKNOWN
from genomics_data_index.storage.model.QueryFeature import QueryFeature
from genomics_data_index.storage.model.QueryFeatureMutation import QueryFeatureMutation


class QueryFeatureMutationSPDI(QueryFeatureMutation):

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
    def id_no_prefix(self) -> Optional[str]:
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

    def has_deletion_sequence(self) -> bool:
        return isinstance(self.ref, str)

    def deletion_length(self) -> int:
        if self.has_deletion_sequence():
            return len(self.deletion)
        else:
            return self.deletion

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
        return self.insertion == NUCLEOTIDE_UNKNOWN

    def to_unknown(self) -> QueryFeature:
        return QueryFeatureMutationSPDI(':'.join([
            self._seq_name,
            str(self._pos),
            str(self._ref),
            NUCLEOTIDE_UNKNOWN
        ]))

    def to_unknown_explode(self) -> List[QueryFeature]:
        if self.deletion_length() == 1:
            unknowns = [self.to_unknown()]
        elif self.has_deletion_sequence():
            unknowns = []
            for i, c in enumerate(self.ref):
                unknowns.append(
                    QueryFeatureMutationSPDI(f'{self.sequence}:{self.position + i}:{c}:{NUCLEOTIDE_UNKNOWN}'))
        else:
            unknowns = []
            for i in range(self.ref):
                unknowns.append(QueryFeatureMutationSPDI(f'{self.sequence}:{self.position + i}:1:{NUCLEOTIDE_UNKNOWN}'))
        return unknowns
