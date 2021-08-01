from __future__ import annotations

from genomics_data_index.storage.model.QueryFeature import QueryFeature
from genomics_data_index.storage.model.QueryFeatureHGVS import QueryFeatureHGVS
from genomics_data_index.storage.model.QueryFeatureHGVSGN import QueryFeatureHGVSGN
from genomics_data_index.storage.model.QueryFeatureMLST import QueryFeatureMLST
from genomics_data_index.storage.model.QueryFeatureMutationSPDI import QueryFeatureMutationSPDI


class QueryFeatureFactory:
    factory_instance = None

    def __init__(self):
        # Maps the prefix string (e.g. 'hgvs:') to a function used to create the specific QueryFeature from the ID
        self._registered_feature_creators = {
            QueryFeatureHGVS.PREFIX: lambda i: QueryFeatureHGVS.create_from_id(i),
            QueryFeatureHGVSGN.PREFIX: lambda i: QueryFeatureHGVSGN.create_from_id(i),
            QueryFeatureMLST.PREFIX: lambda i: QueryFeatureMLST.create_from_id(i),
        }

    def create_feature(self, id: str) -> QueryFeature:
        if id is None:
            raise Exception('id is None')

        for prefix in self._registered_feature_creators:
            if id.startswith(prefix):
                return self._registered_feature_creators[prefix](id)

        # By default assume feature is an SPDI feature
        try:
            return QueryFeatureMutationSPDI(id)
        except Exception as e:
            raise Exception(f'Unknown feature type for feature=[{id}]', e)

    @classmethod
    def instance(cls) -> QueryFeatureFactory:
        if cls.factory_instance is None:
            cls.factory_instance = QueryFeatureFactory()
        return cls.factory_instance
