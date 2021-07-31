import abc
from typing import Dict, Generator, Tuple

from genomics_data_index.api.query.features.FeaturesSummarizer import FeaturesSummarizer
from genomics_data_index.configuration.connector.DataIndexConnection import DataIndexConnection
from genomics_data_index.storage.SampleSet import SampleSet
from genomics_data_index.storage.model.db import FeatureSamples


class FeaturesFromIndexSummarizer(FeaturesSummarizer, abc.ABC):

    def __init__(self, connection: DataIndexConnection):
        super().__init__(connection=connection)

    def _feature_sample_count_iter(self, present_features: Dict[str, FeatureSamples],
                                   present_samples: SampleSet) -> Generator[Tuple[FeatureSamples, int], None, None]:
        for feature_id in present_features:
            feature = present_features[feature_id]
            samples_in_feature = present_samples.intersection(feature.sample_ids)
            yield feature, len(samples_in_feature)
