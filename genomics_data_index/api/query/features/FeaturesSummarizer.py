import abc

import pandas as pd

from genomics_data_index.storage.SampleSet import SampleSet
from genomics_data_index.configuration.connector.DataIndexConnection import DataIndexConnection


class FeaturesSummarizer(abc.ABC):
    FEATURES_SELECTIONS = ['all', 'unique']

    def __init__(self, connection: DataIndexConnection):
        self._connection = connection

    @abc.abstractmethod
    def summary(self, sample_set: SampleSet) -> pd.DataFrame:
        """
        Given a samples query, summarizes the features for all samples in this query.
        :return: A dataframe summarizing features in this set of samples.
        """
        pass
