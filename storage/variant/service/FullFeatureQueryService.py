from typing import List, Any
import abc

import pandas as pd

from storage.variant.service.QueryService import QueryService, QueryFeature
from storage.variant.service import SampleService


class FullFeatureQueryService(QueryService):

    def __init__(self, sample_service: SampleService):
        super().__init__()
        self._sample_service = sample_service

    @abc.abstractmethod
    def get_correct_query_feature(self) -> Any:
        pass

    def validate_query_features(self, features: List[QueryFeature]):
        for feature in features:
            if not isinstance(feature, self.get_correct_query_feature()):
                raise Exception(f'feature=[{feature}] is not of type {self.get_correct_query_feature()}')

    def _find_by_features_internal(self, features: List[QueryFeature], include_unknown: bool) -> pd.DataFrame:
        self.validate_query_features(features)

        feature_samples = self._sample_service.find_samples_by_features(features)

        data = []
        for id in feature_samples:
            for sample in feature_samples[id]:
                data.append([id, sample.name, sample.id, 'Present'])

        results_df = pd.DataFrame(data=data, columns=['Feature', 'Sample Name', 'Sample ID', 'Status'])
        if include_unknown:
            unknown_features_df = self._get_unknown_features(features)
            results_df = pd.concat([results_df, unknown_features_df])

        return results_df.sort_values(['Feature', 'Sample Name'])

    @abc.abstractmethod
    def _get_unknown_features(self, features: List[QueryFeature]) -> pd.DataFrame:
        pass
