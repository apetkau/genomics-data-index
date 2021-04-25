from pathlib import Path
from typing import List, Dict, Any, Set, cast

import pandas as pd

from genomics_data_index.storage.model import MLST_UNKNOWN_ALLELE
from genomics_data_index.storage.model.QueryFeature import QueryFeature
from genomics_data_index.storage.model.QueryFeatureMLST import QueryFeatureMLST
from genomics_data_index.storage.service.FullFeatureQueryService import FullFeatureQueryService
from genomics_data_index.storage.service.MLSTService import MLSTService
from genomics_data_index.storage.service.SampleService import SampleService


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

    def expand_feature(self, feature: QueryFeature) -> List[QueryFeature]:
        if not feature.is_wild():
            return [feature]

        mlst_feature = cast(QueryFeatureMLST, feature)

        new_query_features = []
        if mlst_feature.locus == mlst_feature.WILD:
            loci_with_found_alleles = set()
            all_loci = set()
            all_loci_alleles = self._mlst_service.get_all_loci_alleles(scheme=mlst_feature.scope)
            for (locus, allele) in all_loci_alleles:
                all_loci.add(locus)

                if allele != MLST_UNKNOWN_ALLELE:
                    loci_with_found_alleles.add(locus)
                    new_query_features.append(QueryFeatureMLST.create_feature(scheme=mlst_feature.scope,
                                                                              locus=locus,
                                                                              allele=allele))

            loci_with_only_unknown_alleles = all_loci - loci_with_found_alleles
            # Add feature for loci with invalid alleles if they exist
            for locus in loci_with_only_unknown_alleles:
                new_query_features.append(QueryFeatureMLST.create_feature(scheme=mlst_feature.scope,
                                                                          locus=locus,
                                                                          allele=MLST_UNKNOWN_ALLELE))
        elif mlst_feature.allele == mlst_feature.WILD:
            all_alleles = self._mlst_service.get_all_alleles(scheme=mlst_feature.scope,
                                                             locus=mlst_feature.locus)
            for allele in all_alleles:
                if allele != MLST_UNKNOWN_ALLELE:
                    new_query_features.append(QueryFeatureMLST.create_feature(scheme=mlst_feature.scope,
                                                                              locus=mlst_feature.locus,
                                                                              allele=allele))
        return new_query_features

    def _get_unknown_features(self, features: List[QueryFeature]) -> pd.DataFrame:
        data = []
        unknown_feature_map = {}
        unknown_features_list = []
        for feature in features:
            unknown_feature = feature.to_unknown()
            unknown_features_list.append(unknown_feature)
            unknown_feature_map[unknown_feature.id] = feature

        unknown_feature_samples = self._sample_service.find_samples_by_features(unknown_features_list)
        for uid in unknown_feature_samples:
            for sample in unknown_feature_samples[uid]:
                data.append([unknown_feature_map[uid].id, sample.name, sample.id, 'Unknown'])

        return pd.DataFrame(data, columns=['Feature', 'Sample Name', 'Sample ID', 'Status'])

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
