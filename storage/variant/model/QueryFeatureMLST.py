from storage.variant.model import MLST_UNKNOWN_ALLELE
from storage.variant.model.QueryFeature import QueryFeature


class QueryFeatureMLST(QueryFeature):

    def __init__(self, sla: str):
        super().__init__()
        self._sla = sla

        scheme, locus, allele = self._sla.split(':')
        self._scheme = scheme
        self._locus = locus
        self._allele = allele

    @property
    def id(self):
        return self._sla

    @property
    def scope(self):
        return self._scheme

    @property
    def locus(self):
        return self._locus

    @property
    def allele(self):
        return self._allele

    def to_unknown(self) -> QueryFeature:
        return QueryFeatureMLST(':'.join([
            self._scheme,
            self._locus,
            MLST_UNKNOWN_ALLELE
        ]))