from pathlib import Path
from typing import List, Dict, Any, Set

import pandas as pd

from storage.variant.service.FullFeatureQueryService import FullFeatureQueryService
from storage.variant.service.MLSTService import MLSTService
from storage.variant.service.QueryService import QueryFeature
from storage.variant.service.SampleService import SampleService


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


class MLSTQueryService(FullFeatureQueryService):

    def __init__(self, mlst_service: MLSTService,
                 sample_service: SampleService):
        super().__init__(sample_service=sample_service)
        self._mlst_service = mlst_service

    def get_correct_query_feature(self) -> Any:
        return QueryFeatureMLST

    def _get_feature_scope_sample_counts(self, feature_scopes: Set[str]) -> Dict[str, int]:
        scheme_sample_counts = {
            scope: self._sample_service.count_samples_associated_with_mlst_scheme(scope) for scope in feature_scopes
        }

        return scheme_sample_counts

    def _get_unknown_features(self, features: List[QueryFeature]) -> pd.DataFrame:
        raise Exception('Not implemented')

    def _find_matches_genome_files_internal(self, sample_reads: Dict[str, List[Path]],
                                            distance_threshold: float = None) -> pd.DataFrame:
        raise Exception('Not implemented')

    def _find_matches_internal(self, sample_names: List[str], distance_threshold: float):
        raise Exception('Not implemented')

    def _pairwise_distance_internal(self, samples: List[str]) -> pd.DataFrame:
        raise Exception('Not implemented')

    def _differences_between_genomes_internal(self, sample1: str, sample2: str):
        raise Exception('Not implemented')

    def get_data_type(self) -> str:
        return 'mlst'