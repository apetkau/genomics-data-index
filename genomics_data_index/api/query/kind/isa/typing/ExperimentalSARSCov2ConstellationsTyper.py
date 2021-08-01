import logging

from genomics_data_index.api.query.SamplesQuery import SamplesQuery
from genomics_data_index.api.query.kind.isa.typing.SamplesTypingIsaKind import SamplesTypingIsaKind

logger = logging.getLogger(__name__)


class ExperimentalSARSCov2ConstellationsTyper(SamplesTypingIsaKind):
    """
    This is an experimental/test class to try out typing SARS-CoV-2 based on constellations of variants/mutations.
    Delta mutations derived from https://github.com/cov-lineages/constellations/blob/main/constellations/definitions/cB.1.617.2.json
    In the future I can see about making this more generalized.
    """

    def __init__(self, sequence_name: str = 'MN996528.1'):
        super().__init__()

        self._sequence_name = sequence_name

        initial_mutations = {
            "S:T19R",
            "S:G142D",
            "S:L452R",
            "S:T478K",
            "S:P681R",
            "S:D950N",
            "ORF3a:S26L",
            "M:I82T",
            "ORF7a:V82A",
            "ORF7a:T120I",
            "N:D63G",
            "N:R203M",
            "N:D377Y",
        }
        initial_mutations_split = map(lambda x: x.split(':'), initial_mutations)
        self._mutation_ids = {f'hgvs_gn:{self._sequence_name}:{m[0]}:p.{m[1]}' for m in initial_mutations_split}
        logger.info(f'mutation_ids={self._mutation_ids}')

    @property
    def name(self) -> str:
        return 'constellations_typer'

    @property
    def version(self) -> str:
        return 'experimental'

    def isa_type(self, data: str, query: SamplesQuery) -> SamplesQuery:
        # Initial implementation, I just look for what has all these mutations
        q = query
        for mutation in self._mutation_ids:
            q = q.hasa(mutation, kind='mutation')

        return q
