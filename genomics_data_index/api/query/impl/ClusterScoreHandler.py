import pandas as pd

from genomics_data_index.api.query.impl.TreeSamplesQuery import TreeSamplesQuery
from genomics_data_index.configuration.connector.DataIndexConnection import DataIndexConnection
from genomics_data_index.storage.SampleSet import SampleSet


class ClusterScoreHandler:

    def __init__(self, samples_query: TreeSamplesQuery):
        self._samples_query = samples_query

    def score_clusters(self, cluster_series: pd.Series) -> pd.Series:
        cluster_sample_sets = cluster_series.groupby().agg({'Sample IDs': SampleSet})
