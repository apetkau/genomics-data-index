from __future__ import annotations

from typing import Optional

from genomics_data_index.storage.model.QueryFeatureMutation import QueryFeatureMutation
from genomics_data_index.storage.model.QueryFeature import QueryFeature


class QueryFeatureHGVS(QueryFeatureMutation):

    def __init__(self, sequence_name: str, gene: Optional[str], variant: str):
        super().__init__()

        self._sequence_name = sequence_name
        self._gene = gene
        self._variant = variant

        if self._variant is not None:
            if self._gene is None:
                self._hgvs_id = f'{sequence_name}:{variant}'
            else:
                self._hgvs_id = f'{sequence_name}:{gene}:{variant}'
        else:
            self._hgvs_id = None

    def has_gene(self) -> bool:
        return self._gene is not None

    def has_id(self) -> bool:
        return self._hgvs_id is not None

    def is_nucleotide(self) -> bool:
        return self.has_id() and (self._variant.startswith('n.') or self._variant.startswith('c.'))

    def is_protein(self) -> bool:
        return self.has_id() and self._variant.startswith('p.')

    @property
    def sequence(self) -> str:
        return self._sequence_name

    @property
    def gene(self) -> str:
        return self._gene

    @property
    def variant(self) -> str:
        return self._variant

    @property
    def id(self) -> Optional[str]:
        if self.has_id():
            return f'hgvs:{self._hgvs_id}'
        else:
            return None

    def is_unknown(self) -> bool:
        raise NotImplementedError('Not implemented')

    @property
    def scope(self) -> str:
        return self.sequence

    def to_unknown(self) -> QueryFeature:
        raise NotImplementedError('Not implemented')

    @classmethod
    def create(cls, sequence_name: str, gene: Optional[str], variant: str) -> QueryFeatureHGVS:
        # If variant is given as nucleotide coordinates with respect to the sequence 'n.' then ignore gene id
        if variant is not None and variant.startswith('n.'):
            gene = None

        return QueryFeatureHGVS(sequence_name=sequence_name, gene=gene, variant=variant)

    @classmethod
    def create_from_id(cls, hgvs_id: str) -> QueryFeatureHGVS:
        """
        Creates a new HGVS feature. The identifier is in the format hgvs:[sequence]:[gene]:[variant] or
        hgvs:[sequence]:[variant]. This roughly corresponds to the HGVS format <http://varnomen.hgvs.org/>.
        Specifically, the [gene]:[variant] (or [reference]:[variant]) part corresponds to the HGVS identifier.
        I prefix the identifier with the string 'hgvs:' to indicate it is an HGVS identifier. I also include
        the reference genome name in addition to the gene name in cases where the variant is given in gene coordinates.
        :param hgvs_id: The (modified) HGVS identifier.
        """
        if hgvs_id.startswith(f'hgvs{cls.SPLIT_CHAR}'):
            hgvs_id_strip = hgvs_id[len(f'hgvs{cls.SPLIT_CHAR}'):]
        else:
            hgvs_id_strip = hgvs_id

        tokens = hgvs_id_strip.split(cls.SPLIT_CHAR)
        if len(tokens) == 2:
            sequence_name = tokens[0]
            gene = None
            variant = tokens[1]
        elif len(tokens) == 3:
            sequence_name = tokens[0]
            gene = tokens[1]
            variant = tokens[2]
        else:
            raise Exception(f'Invalid number of items in hgvs_id=[{hgvs_id}].'
                            f' Should be in the form [hgvs:reference:gene:variant] or [hgvs:reference:variant].')

        return QueryFeatureHGVS(sequence_name=sequence_name, gene=gene, variant=variant)
