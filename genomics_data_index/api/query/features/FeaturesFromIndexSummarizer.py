import abc
from typing import Dict, Generator, Tuple, List, Any

from genomics_data_index.api.query.features.FeaturesSummarizer import FeaturesSummarizer
from genomics_data_index.configuration.connector.DataIndexConnection import DataIndexConnection
from genomics_data_index.storage.SampleSet import SampleSet
from genomics_data_index.storage.model.db import FeatureSamples


class FeaturesFromIndexSummarizer(FeaturesSummarizer, abc.ABC):

    def __init__(self, connection: DataIndexConnection):
        super().__init__(connection=connection)

    def _create_feature_sample_count_rows(self, present_features: Dict[str, FeatureSamples],
                                   present_samples: SampleSet) -> List[Any]:
        data = []
        total = len(present_samples)
        for feature_id in present_features:
            feature = present_features[feature_id]
            samples_in_feature = present_samples.intersection(feature.sample_ids)
            sample_count = len(samples_in_feature)
            if sample_count > 0:
                data.append(self._create_feature_sample_count_row(feature,
                                                                  sample_count=sample_count,
                                                                  total=total))
        return data

    @abc.abstractmethod
    def _create_feature_sample_count_row(self, feature: FeatureSamples,
                                         sample_count: int, total: int) -> List[Any]:
        pass
