from __future__ import annotations
from typing import List

import pandas as pd

from storage.connector.DataIndexConnection import DataIndexConnection


class GenomicDataStore:

    def __init__(self, connection: DataIndexConnection):
        self._connection = connection

    def count_samples(self) -> int:
        return self._connection.sample_service.count_samples()

    def sample_names(self) -> List[str]:
        return [s.name for s in self._connection.sample_service.get_samples()]

    def count_references(self) -> int:
        return self._connection.reference_service.count_reference_genomes()

    def reference_names(self) -> List[str]:
        return [r.name for r in self._connection.reference_service.get_reference_genomes()]

    def count_mutations(self, reference_genome: str, include_unknown: bool = False) -> int:
        return self._connection.variation_service.count_on_reference(reference_genome,
                                                                     include_unknown=include_unknown)

    def mutations_summary(self, reference_genome: str, include_unknown: bool = False) -> pd.Series:
        mutation_counts = self._connection.variation_service.mutation_counts_on_reference(reference_genome,
                                                                                          include_unknown=include_unknown)

        return pd.Series(mutation_counts, name='Mutations')
