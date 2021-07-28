from __future__ import annotations

from typing import List, Optional

from genomics_data_index.storage.model.QueryFeature import QueryFeature
from genomics_data_index.storage.model.QueryFeatureHGVS import QueryFeatureHGVS


class QueryFeatureHGVSGN(QueryFeatureHGVS):
    PREFIX = 'hgvs_gn:'

    def __init__(self, sequence_name: str, gene: str, variant: str):
        super().__init__(sequence_name=sequence_name, gene=gene, variant=variant)

    def is_unknown(self) -> bool:
        raise NotImplementedError('Not implemented')

    def to_unknown(self) -> QueryFeature:
        raise NotImplementedError('Not implemented')

    def to_unknown_explode(self) -> List[QueryFeature]:
        raise NotImplementedError('Not implemented')

    @classmethod
    def create(cls, sequence_name: str, gene: Optional[str], variant: Optional[str]) -> QueryFeatureHGVSGN:
        # If variant is given as nucleotide coordinates with respect to the sequence 'n.' then ignore gene id
        if variant is not None and variant.startswith('n.'):
            gene = None

        return QueryFeatureHGVSGN(sequence_name=sequence_name, gene=gene, variant=variant)

    @classmethod
    def create_from_id(cls, hgvs_id: str, prefix: str = None) -> QueryFeatureHGVSGN:
        """
        Creates a new HGVSGN feature which uses the gene name (instead of locus id) for the "gene" component. This
        is nearly identical to an HGVS identifier, except that the [gene] portion in the id
        (hgvs_gn:[sequence]:[gene]:[variant]) now corresponds to the gene name instead of locus identifier. For example,
        hgvs_gn:MN996528.1:S:p.D614G to specify the D614G mutation on the spike (S) protein name instead of the
        corresponding locus identifier.
        :param hgvs_id: The (modified) HGVS identifier.
        :param prefix: The prefix to use for this identifier. Used to override prefix for HGVSGN
        :return: A new HGVS query feature.
        """
        if prefix is None:
            prefix = cls.PREFIX

        sequence_name, gene, variant = QueryFeatureHGVS.split_id(hgvs_id, prefix=prefix)
        return QueryFeatureHGVSGN(sequence_name=sequence_name, gene=gene, variant=variant)
