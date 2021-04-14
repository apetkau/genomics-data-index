from __future__ import annotations

from storage.variant.model import MLST_UNKNOWN_ALLELE
from storage.variant.model.QueryFeature import QueryFeature


class QueryFeatureMLST(QueryFeature):

    def __init__(self, sla: str):
        super().__init__()
        self._sla = sla

        scheme, locus, allele = self._sla.split(self.SPLIT_CHAR)

        if scheme == self.WILD:
            raise Exception(f'Cannot set seq to be wild ({self.WILD}): {sla}')
        else:
            self._scheme = scheme

        if locus == self.WILD:
            raise Exception(f'Cannot set seq to be wild ({self.WILD}): {sla}')
        else:
            self._locus = locus

        if allele is None:
            self._allele = self.WILD
        else:
            self._allele = allele

    @property
    def id(self) -> str:
        return self._sla

    @property
    def scope(self) -> str:
        return self._scheme

    @property
    def locus(self) -> str:
        return self._locus

    @property
    def allele(self) -> str:
        return self._allele

    @classmethod
    def create_feature(cls, scheme: str, locus: str, allele: str) -> QueryFeatureMLST:
        return QueryFeatureMLST(cls.SPLIT_CHAR.join([scheme, locus, allele]))

    def to_unknown(self) -> QueryFeature:
        return QueryFeatureMLST(':'.join([
            self._scheme,
            self._locus,
            MLST_UNKNOWN_ALLELE
        ]))