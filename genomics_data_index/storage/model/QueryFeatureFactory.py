from __future__ import annotations

from genomics_data_index.storage.model.QueryFeature import QueryFeature
from genomics_data_index.storage.model.QueryFeatureHGVS import QueryFeatureHGVS


class QueryFeatureFactory:
    factory_instance = None

    registered_feature_creators = {
        QueryFeatureHGVS.PREFIX: lambda i: QueryFeatureHGVS.create_from_id(i),
    }

    def __init__(self):
        pass

    def create_feature(self, id: str) -> QueryFeature:
        if id is None:
            raise Exception('id is None')

        for prefix in self.registered_feature_creators:
            if id.startswith(prefix):
                return self.registered_feature_creators[prefix](id)

        raise Exception(f'Unknown feature type for feature=[{id}]')

    @classmethod
    def instance(cls) -> QueryFeatureFactory:
        if cls.factory_instance is None:
            cls.factory_instance = QueryFeatureFactory()
        return cls.factory_instance
