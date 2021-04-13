from storage.variant.model import NUCLEOTIDE_UNKNOWN
from storage.variant.model.QueryFeature import QueryFeature


class QueryFeatureMutation(QueryFeature):

    def __init__(self, spdi: str):
        super().__init__()
        self._spdi = spdi

        seq, pos, ref, alt = self._spdi.split(':')
        self._seq_name = seq
        self._pos = int(pos)
        self._ref = ref
        self._alt = alt

    @property
    def id(self):
        return self._spdi

    @property
    def scope(self):
        return self._seq_name

    @property
    def position(self):
        return self._pos

    @property
    def start(self):
        return self.position

    @property
    def stop(self):
        return self.position + len(self.ref)

    @property
    def ref(self):
        return self._ref

    @property
    def alt(self):
        return self._alt

    def to_unknown(self) -> QueryFeature:
        return QueryFeatureMutation(':'.join([
            self._seq_name,
            self._pos,
            self._ref,
            NUCLEOTIDE_UNKNOWN
        ]))