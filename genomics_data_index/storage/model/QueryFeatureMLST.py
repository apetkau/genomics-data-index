from __future__ import annotations

from typing import List, Optional

from genomics_data_index.storage.model import MLST_UNKNOWN_ALLELE
from genomics_data_index.storage.model.QueryFeature import QueryFeature


class QueryFeatureMLST(QueryFeature):
    PREFIX = 'mlst:'

    def __init__(self, mlst_id: str):
        super().__init__()

        if mlst_id.startswith(QueryFeatureMLST.PREFIX):
            sla = mlst_id[len(QueryFeatureMLST.PREFIX):]
        else:
            raise Exception(f'mlst_id=[{mlst_id}] must start with [{QueryFeatureMLST.PREFIX}]')

        values = sla.split(self.SPLIT_CHAR)
        if len(values) == 3:
            scheme, locus, allele = values[0], values[1], values[2]
        elif len(values) == 2:
            scheme, locus = values[0], values[1]
            allele = None
        else:
            raise Exception(
                f'Invalid MLST feature [{sla}]. Must be either mlst:scheme:locus or mlst:scheme:locus:allele')

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
            raise Exception(
                f'Unsupported to set wild ({self.WILD}) for locus and not for allele ({self.WILD}): {mlst_id}')

    @property
    def id_no_prefix(self) -> Optional[str]:
        return f'{self.scheme}:{self.locus}:{self.allele}'

    @property
    def prefix(self) -> str:
        return self.PREFIX

    @property
    def scope(self) -> str:
        return self.scheme

    @property
    def scheme(self) -> str:
        return self._scheme

    @property
    def locus(self) -> str:
        return self._locus

    @property
    def allele(self) -> str:
        return self._allele

    def is_unknown(self) -> bool:
        return self.allele == MLST_UNKNOWN_ALLELE

    @classmethod
    def create_feature(cls, scheme: str, locus: str, allele: str) -> QueryFeatureMLST:
        return QueryFeatureMLST(cls.SPLIT_CHAR.join([scheme, locus, allele]))

    @classmethod
    def to_query_id(cls, mlst_id: str) -> str:
        if not mlst_id.startswith(QueryFeatureMLST.PREFIX):
            mlst_id = 'mlst:' + mlst_id
        return QueryFeatureMLST.create_from_id(mlst_id).id

    @classmethod
    def create_from_id(cls, mlst_id: str) -> QueryFeatureMLST:
        return QueryFeatureMLST(mlst_id)

    def to_unknown(self) -> QueryFeature:
        return QueryFeatureMLST(f'{self.prefix}{self.scheme}:{self.locus}:{MLST_UNKNOWN_ALLELE}')

    def to_unknown_explode(self) -> List[QueryFeature]:
        return [self.to_unknown()]
