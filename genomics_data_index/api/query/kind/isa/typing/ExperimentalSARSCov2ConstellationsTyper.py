import logging
from collections import Counter

from genomics_data_index.api.query.SamplesQuery import SamplesQuery
from genomics_data_index.api.query.kind.isa.typing.SamplesTypingIsaKind import SamplesTypingIsaKind
from genomics_data_index.storage.SampleSet import SampleSet

logger = logging.getLogger(__name__)


class ExperimentalSARSCov2ConstellationsTyper(SamplesTypingIsaKind):
    """
    This is an experimental/test class to try out typing SARS-CoV-2 based on constellations of variants/mutations.
    Delta mutations derived from https://github.com/cov-lineages/constellations/blob/main/constellations/definitions/cB.1.617.2.json
    In the future I can see about making this more generalized.
    """

    def __init__(self, sequence_name: str = 'MN996528.1', min_alt: int = 5):
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

        self._number_mutations = len(self._mutation_ids)
        self._min_alt = min_alt
        self._max_ref = 3

    @property
    def name(self) -> str:
        return 'constellations_typer'

    @property
    def version(self) -> str:
        return 'experimental'

    def perfect_matches(self, data: str, query: SamplesQuery) -> SamplesQuery:
        q = query
        for mutation in self._mutation_ids:
            q = q.hasa(mutation, kind='mutation')

        return q

    def imperfect_matches(self, data: str, query: SamplesQuery) -> SamplesQuery:
        # Use a Counter to keep track of number of mutations associated with each sample
        samples_present_counter = Counter(query.sample_set)
        # sample_unknown_mutation_counter = Counter(query.unknown_set)

        for mutation in self._mutation_ids:
            mutation_query = query.hasa(mutation, kind='mutation')
            samples_present_counter.update(mutation_query.sample_set)
            # sample_unknown_mutation_counter.update(mutation_query.unknown_set)

        sample_ids_pass = {s for s in samples_present_counter if samples_present_counter[s] > self._min_alt}
        sample_set_present = SampleSet(sample_ids_pass)

        return query.intersect(sample_set_present)

    def isa_type(self, data: str, query: SamplesQuery) -> SamplesQuery:
        # I first look for perfect matches to all mutations to avoid having to count mutations for these
        query_perfect_match = self.perfect_matches(data, query)
        query_imperfect_matches = query.reset_universe() & (~(query_perfect_match.reset_universe()))

        # Now I have to look for imperfect matches to mutations and include them in my final result
        query_imperfect_matches = self.imperfect_matches(data, query_imperfect_matches)

        return query_perfect_match | query_imperfect_matches
