from typing import Generator, Tuple, Dict, cast

import pandas as pd

from genomics_data_index.api.query.features.FeaturesSummarizer import FeaturesSummarizer
from genomics_data_index.configuration.connector.DataIndexConnection import DataIndexConnection
from genomics_data_index.storage.service.MLSTService import MLSTService
from genomics_data_index.storage.SampleSet import SampleSet
from genomics_data_index.storage.model.db import FeatureSamples, MLSTAllelesSamples


class MLSTFeaturesSummarizer(FeaturesSummarizer):

    def __init__(self, connection: DataIndexConnection, scheme: str = None):
        super().__init__(connection=connection)
        self._scheme = scheme

    def _feature_sample_count_iter(self, present_features: Dict[str, FeatureSamples],
                                   present_samples: SampleSet) -> Generator[Tuple[FeatureSamples, int], None, None]:
        for feature_id in present_features:
            feature = present_features[feature_id]
            samples_in_feature = present_samples.intersection(feature.sample_ids)
            yield feature, len(samples_in_feature)

    def summary(self, present_samples: SampleSet,
                unknown_samples: SampleSet,
                absent_samples: SampleSet,
                selection: str) -> pd.DataFrame:
        mlst_service: MLSTService = self._connection.mlst_service

        # TODO: This isn't the best implementation right now since I have to load
        # all MLST feature => sample mappings and then iterate over them
        # To do this more efficiently I would probably need to modify my data model
        mlst_present_features = mlst_service.get_features(scheme=self._scheme, include_unknown=False)
        data = []
        for feature, sample_count in self._feature_sample_count_iter(present_features=mlst_present_features,
                                                                     present_samples=present_samples):
            if sample_count > 0:
                feature = cast(MLSTAllelesSamples, feature)
                data.append([feature.query_id, feature.scheme, feature.locus, feature.allele, sample_count])
        summary_df = pd.DataFrame(data,
                                  columns=['MLST Feature', 'Scheme', 'Locus', 'Allele', 'Count'])
        summary_df['Total'] = len(present_samples)
        summary_df['Percent'] = 100 * (summary_df['Count'] / summary_df['Total'])

        return summary_df.set_index('MLST Feature')
