from pathlib import Path
from typing import List, Dict

import pandas as pd

from storage.variant.service.QueryService import QueryService


class KmerQueryService(QueryService):

    def __init__(self):
        super().__init__()

    def _find_matches_genome_files_internal(self, sample_reads: Dict[str, List[Path]],
                                            distance_threshold: float = None) -> pd.DataFrame:
        raise Exception('Not implemented')

    def _find_matches_internal(self, samples, distance_threshold: float):
        raise Exception('Not implemented')

    def _pairwise_distance_internal(self, samples: List[str]) -> pd.DataFrame:
        raise Exception('Not implemented')

    def _differences_between_genomes_internal(self, sample1: str, sample2: str):
        raise Exception('Not implemented')

    def get_data_type(self) -> str:
        return 'kmer'
