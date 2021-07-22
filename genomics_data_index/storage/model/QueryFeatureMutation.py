import abc

from genomics_data_index.storage.model.QueryFeature import QueryFeature


class QueryFeatureMutation(QueryFeature, abc.ABC):

    def __init__(self):
        super().__init__()
