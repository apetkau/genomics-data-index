from pathlib import Path
from typing import List, Dict

import pandas as pd

from storage.variant.service.QueryService import QueryService, QueryFeature
from storage.variant.service.SampleService import SampleService
from storage.variant.service.MLSTService import MLSTService


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
    def scheme(self):
        return self._scheme

    @property
    def locus(self):
        return self._locus

    @property
    def allele(self):
        return self._allele


class MLSTQueryService(QueryService):

    def __init__(self, mlst_service: MLSTService,
                 sample_service: SampleService):
        super().__init__()
        self._mlst_service = mlst_service
        self._sample_service = sample_service

    def _find_matches_genome_files_internal(self, sample_reads: Dict[str, List[Path]],
                                            distance_threshold: float = None) -> pd.DataFrame:
        pass

    def _find_matches_internal(self, sample_names: List[str], distance_threshold: float):
        pass

    def _pairwise_distance_internal(self, samples: List[str]) -> pd.DataFrame:
        pass

    def _differences_between_genomes_internal(self, sample1: str, sample2: str):
        pass

    def get_data_type(self) -> str:
        return 'mlst'
