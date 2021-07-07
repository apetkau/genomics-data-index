from abc import ABC

from genomics_data_index.storage.model.QueryFeature import QueryFeature


class QueryFeatureMutation(QueryFeature, ABC):

    def __init__(self):
        super().__init__()
