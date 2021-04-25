from __future__ import annotations

from genomics_data_index.variant.model import MLST_UNKNOWN_ALLELE
from genomics_data_index.variant.model.QueryFeature import QueryFeature


class QueryFeatureMLST(QueryFeature):

    def __init__(self, sla: str):
        super().__init__()
        self._sla = sla

        scheme = None
        locus = None
        allele = None

        values = self._sla.split(self.SPLIT_CHAR)
        if len(values) == 3:
            scheme, locus, allele = values[0], values[1], values[2]
        elif len(values) == 2:
            scheme, locus = values[0], values[1]
        else:
            raise Exception(f'Invalid MLST feature [{sla}]. Must be either scheme:locus or scheme:locus:allele')

        if scheme == self.WILD:
            raise Exception(f'Cannot set seq to be wild ({self.WILD}): {sla}')
        else:
            self._scheme = scheme

        if locus is None:
            self._locus = self.WILD
        else:
            self._locus = locus

        if allele is None:
            self._allele = self.WILD
        else:
            self._allele = allele

        if self._locus == self.WILD and self._allele != self.WILD:
            raise Exception(f'Unsupported to set wild ({self.WILD}) for locus and not for allele ({self.WILD}): {sla}')

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

    def is_unknown(self) -> bool:
        return self._allele == MLST_UNKNOWN_ALLELE

    @classmethod
    def create_feature(cls, scheme: str, locus: str, allele: str) -> QueryFeatureMLST:
        return QueryFeatureMLST(cls.SPLIT_CHAR.join([scheme, locus, allele]))

    def to_unknown(self) -> QueryFeature:
        return QueryFeatureMLST(':'.join([
            self._scheme,
            self._locus,
            MLST_UNKNOWN_ALLELE
        ]))
