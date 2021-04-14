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
