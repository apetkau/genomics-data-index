import abc

import pandas as pd

from genomics_data_index.storage.SampleSet import SampleSet


class FeaturesSummarizer(abc.ABC):
    FEATURES_SELECTIONS = ['all', 'unique']

    def __init__(self):
        pass

    @abc.abstractmethod
    def summary(self, present_samples: SampleSet,
                unknown_samples: SampleSet,
                absent_samples: SampleSet,
                selection: str) -> pd.DataFrame:
        """
        Given a samples query, summarizes the features for all samples in this query.
        :param present_samples: The samples that are present.
        :param unknown_samples: The samples that are unknown.
        :param absent_samples: The samples not present (or unknown).
        :param selection: THe selection of samples to summarize ('all' or 'unique').
        :return: A dataframe summarizing features in this query.
        """
        pass
