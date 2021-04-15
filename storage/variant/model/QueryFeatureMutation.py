from __future__ import annotations
from typing import Union, Tuple
import re

from storage.variant.model import NUCLEOTIDE_UNKNOWN
from storage.variant.model.QueryFeature import QueryFeature


class QueryFeatureMutation(QueryFeature):

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
    def position(self) -> int:
        return self._pos

    @property
    def start(self) -> int:
        return self.position

    @property
    def stop(self) -> int:
        return self.position + len(self.ref)

    @property
    def ref(self) -> str:
        return self._ref

    @property
    def alt(self) -> str:
        return self._alt

    def is_unknown(self) -> bool:
        return False

    def to_unknown(self) -> QueryFeature:
        return QueryFeatureMutation(':'.join([
            self._seq_name,
            self._pos,
            self._ref,
            NUCLEOTIDE_UNKNOWN
        ]))

    @classmethod
    def convert_deletion(cls, deletion: Union[str, int]) -> int:
        if isinstance(deletion, str):
            if re.match('^\d+$', deletion):
               deletion = int(deletion)
            else:
                if not set(deletion).issubset({'A','T','C','G'}):
                    raise Exception('Deletion must either be an integer or a string with alphabet {A,T,C,G}'
                                    f': {deletion}')
                deletion = len(deletion)
        elif isinstance(deletion, int):
            if deletion < 0:
                raise Exception(f'ref=[{deletion}] must be a non-negative integer')
        else:
            raise Exception(f'ref=[{deletion}] must be either a string or a non-negative integer')

        return deletion

    @classmethod
    def from_spdi(cls, spdi: str) -> Tuple[str, int, int, str]:
        if spdi is None:
            raise Exception('Cannot parse value spdi=None')

        values = spdi.split(':')
        if len(values) != 4:
            raise Exception(f'Incorrect number of items for spdi=[{spdi}]')
        else:
            position = int(values[1])
            deletion = cls.convert_deletion(values[2])

            if position < 0:
                raise Exception(f'Position must be non-negative: {position}')

            return str(values[0]), position, deletion, str(values[3])

    @classmethod
    def to_spdi(cls, sequence_name: str, position: int, ref: Union[str, int], alt: str) -> str:
        if position < 0:
            raise Exception(f'Position must be non-negative: {position}')

        ref = cls.convert_deletion(ref)
        return f'{sequence_name}:{position}:{ref}:{alt}'

    def standardize_feature(self) -> QueryFeatureMutation:
        standardized_spdi = QueryFeatureMutation.to_spdi(self.scope, self.position, self.ref, self.alt)
        return QueryFeatureMutation(standardized_spdi)
