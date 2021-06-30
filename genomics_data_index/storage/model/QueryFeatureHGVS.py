from __future__ import annotations

from typing import Union

from genomics_data_index.storage.model.QueryFeature import QueryFeature


class QueryFeatureHGVS(QueryFeature):

    def __init__(self, hgvs_id: str):
        super().__init__()

        if hgvs_id.startswith(f'hgvs{self.SPLIT_CHAR}'):
            hgvs_id_strip = hgvs_id[len(f'hgvs{self.SPLIT_CHAR}'):]
        else:
            hgvs_id_strip = hgvs_id

        self._hgvs_id = hgvs_id_strip

        tokens = self._hgvs_id.split(self.SPLIT_CHAR)
        if len(tokens) == 2:
            sequence = tokens[0]
            gene = None
            mutation = tokens[1]
        elif len(tokens) == 3:
            sequence = tokens[0]
            gene = tokens[1]
            mutation = tokens[2]
        else:
            raise Exception(f'Invalid number of items in hgvs_id=[{hgvs_id}].'
                            f' Should be in the form [hgvs:sequence:gene:mutation] or [hgvs:sequence:mutation].')

        self._sequence = sequence
        self._gene = gene
        self._mutation = mutation

    def has_gene(self) -> bool:
        return self._gene is not None

    def is_nucleotide(self) -> bool:
        return self._mutation.startswith('n.') or self._mutation.startswith('c.')

    def is_protein(self) -> bool:
        return self._mutation.startswith('p.')

    @property
    def sequence(self) -> str:
        return self._sequence

    @property
    def gene(self) -> str:
        return self._gene

    @property
    def mutation(self) -> str:
        return self._mutation

    @property
    def id(self) -> str:
        return f'hgvs:{self._hgvs_id}'

    def is_unknown(self) -> bool:
        raise NotImplementedError('Not implemented')

    @property
    def scope(self) -> str:
        return self.sequence

    def to_unknown(self) -> QueryFeature:
        raise NotImplementedError('Not implemented')
