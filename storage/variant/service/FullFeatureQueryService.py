import abc
from typing import List, Any, Dict, Set

import pandas as pd

from storage.variant.service import SampleService
from storage.variant.service.QueryService import QueryService
from storage.variant.model.QueryFeature import QueryFeature


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

    def expand_features(self, features: List[QueryFeature]) -> List[QueryFeature]:
        expanded_features = []
        for feature in features:
            if feature.is_wild():
                expanded = self.expand_feature(feature)
                expanded_features.extend(expanded)
            else:
                expanded_features.append(feature)

        return expanded_features

    @abc.abstractmethod
    def expand_feature(self, feature: QueryFeature) -> List[QueryFeature]:
        pass

    def _find_by_features_internal(self, features: List[QueryFeature], include_unknown: bool) -> pd.DataFrame:
        self.validate_query_features(features)
        features = self.expand_features(features)

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
    def _get_feature_scope_sample_counts(self, feature_scopes: Set[str]) -> Dict[str, int]:
        pass

    def _count_by_features_internal(self, features: List[QueryFeature], include_unknown: bool) -> pd.DataFrame:
        self.validate_query_features(features)

        feature_sample_counts = self._sample_service.count_samples_by_features(features)

        scope_set = {f.scope for f in features}
        feature_scope_sample_counts = self._get_feature_scope_sample_counts(scope_set)

        if include_unknown:
            unknown_features_df = self._get_unknown_features(features)
            grouped_df = unknown_features_df.groupby('Feature').agg('count')
            unknown_counts = {}
            absent_counts = {}
            for feature in features:
                if feature.id in grouped_df.index:
                    unknown_counts[feature.id] = grouped_df.loc[feature.id].values[0]
                    absent_counts[feature.id] = feature_scope_sample_counts[feature.scope] \
                                                - feature_sample_counts[feature.id] - unknown_counts[feature.id]
                else:
                    unknown_counts[feature.id] = 0
                    absent_counts[feature.id] = feature_scope_sample_counts[feature.scope] - feature_sample_counts[
                        feature.id]
        else:
            unknown_counts = {f.id: pd.NA for f in features}
            absent_counts = {f.id: feature_scope_sample_counts[f.scope] - feature_sample_counts[f.id]
                             for f in features}

        data = [{
            'Feature': f.id,
            'Present': feature_sample_counts[f.id],
            'Absent': absent_counts[f.id],
            'Unknown': unknown_counts[f.id],
            'Total': feature_scope_sample_counts[f.scope]
        } for f in features]

        count_df = pd.DataFrame(data=data)
        count_df['% Present'] = 100 * count_df['Present'] / count_df['Total']
        count_df['% Absent'] = 100 * count_df['Absent'] / count_df['Total']
        count_df['% Unknown'] = 100 * count_df['Unknown'] / count_df['Total']

        return count_df.sort_values('Feature')

    @abc.abstractmethod
    def _get_unknown_features(self, features: List[QueryFeature]) -> pd.DataFrame:
        pass
