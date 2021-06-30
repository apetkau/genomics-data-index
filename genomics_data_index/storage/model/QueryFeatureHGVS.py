from __future__ import annotations

from typing import Optional

from genomics_data_index.storage.model.QueryFeature import QueryFeature


class QueryFeatureHGVS(QueryFeature):

    def __init__(self, hgvs_id: str):
        """
        Creates a new HGVS feature. The identifier is in the format hgvs:[sequence]:[gene]:[variant] or
        hgvs:[sequence]:[variant]. This roughly corresponds to the HGVS format <http://varnomen.hgvs.org/>.
        Specifically, the [gene]:[variant] (or [reference]:[variant]) part corresponds to the HGVS identifier.
        I prefix the identifier with the string 'hgvs:' to indicate it is an HGVS identifier. I also include
        the reference genome name in addition to the gene name in cases where the variant is given in gene coordinates.
        :param hgvs_id: The (modified) HGVS identifier.
        """

        super().__init__()

        if hgvs_id.startswith(f'hgvs{self.SPLIT_CHAR}'):
            hgvs_id_strip = hgvs_id[len(f'hgvs{self.SPLIT_CHAR}'):]
        else:
            hgvs_id_strip = hgvs_id

        self._hgvs_id = hgvs_id_strip

        tokens = self._hgvs_id.split(self.SPLIT_CHAR)
        if len(tokens) == 2:
            reference_name = tokens[0]
            gene = None
            variant = tokens[1]
        elif len(tokens) == 3:
            reference_name = tokens[0]
            gene = tokens[1]
            variant = tokens[2]
        else:
            raise Exception(f'Invalid number of items in hgvs_id=[{hgvs_id}].'
                            f' Should be in the form [hgvs:reference:gene:variant] or [hgvs:reference:variant].')

        self._reference_name = reference_name
        self._gene = gene
        self._variant = variant

    def has_gene(self) -> bool:
        return self._gene is not None

    def is_nucleotide(self) -> bool:
        return self._variant.startswith('n.') or self._variant.startswith('c.')

    def is_protein(self) -> bool:
        return self._variant.startswith('p.')

    @property
    def reference(self) -> str:
        return self._reference_name

    @property
    def gene(self) -> str:
        return self._gene

    @property
    def variant(self) -> str:
        return self._variant

    @property
    def id(self) -> str:
        return f'hgvs:{self._hgvs_id}'

    def is_unknown(self) -> bool:
        raise NotImplementedError('Not implemented')

    @property
    def scope(self) -> str:
        return self.reference

    def to_unknown(self) -> QueryFeature:
        raise NotImplementedError('Not implemented')

    @classmethod
    def create(cls, sequence_name: str, gene_name: Optional[str], variant: str):
        if gene_name is None:
            hgvs_id = f'hgvs:{sequence_name}:{variant}'
        else:
            hgvs_id = f'hgvs:{sequence_name}:{gene_name}:{variant}'

        return QueryFeatureHGVS(hgvs_id=hgvs_id)
